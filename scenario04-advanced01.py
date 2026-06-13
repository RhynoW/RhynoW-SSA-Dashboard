#!/usr/bin/env python3
"""
scenario04-advanced01.py — 太空態勢儀表板（進階版）
====================================================
Port 5011

以 scenario04-3.py 為基礎，新增三項進階功能：

1-A  向量化 SGP4（SatrecArray）
     - 使用 sgp4.api.SatrecArray 對整個過濾群組批次傳播
     - 一次呼叫 C 層，省去 Python 迴圈呼叫開銷
     - Starlink 全量傳播從 ~10 s 降至 <1 s
     - 自動 fallback：若 SatrecArray 不可用則退回逐顆傳播

2-A  即時近距離配對掃描（KD-tree Proximity Scan）
     - /api/conjunctions：批次傳播全目錄→ ECI 空間建 cKDTree
       → query_pairs 找出當下位置距離 ≤ threshold_km 的配對
     - 輕量版 SSA 指標（非完整 Stage A/B/C Pc 計算）
     - 頂欄顯示「近距離配對 N」計數，點選進入配對清單
     - 點選配對 → 在地球上標出兩星並繪連線，flyTo 中點
     - 詳細接近事件分析請使用 conjunction_pipeline.py

4-A  物件搜尋框
     - 面板頂部輸入框：支援 NORAD ID 精確比對或衛星名稱子串搜尋
     - 即時 debounce（400 ms）→ /api/search → 顯示 top-20
     - 點選結果：在地球上標出該衛星，flyTo，開啟 infoBox

前置準備（與 scenario04-3.py 相同）：
  1. data/globe_texture.jpg  — NASA Blue Marble
  2. data/cesium/            — CesiumJS 1.114 Build/Cesium/
  3. data/borders.geojson    — Natural Earth 國界

額外相依：
  - scipy（KD-tree，已在 conjunction_pipeline.py 使用）
  - sgp4 >= 2.0（SatrecArray 支援）

API：
  GET /api/stats                       → 統計
  GET /api/positions                   → ?ftype=...&fval=...（向量化）
  GET /api/position/<norad_id>         → 單顆衛星位置（供搜尋用）
  GET /api/conjunctions                → ?threshold_km=10&max_pairs=200
  GET /api/search                      → ?q=<NORAD ID 或名稱>
  GET /api/globe_texture               → 離線地球貼圖
  GET /api/layers/borders              → 全球國界 GeoJSON
  GET /api/layers/ssn_stations         → SSN 地面觀測站
  GET /cesium/<path>                   → 本機 CesiumJS
"""

from __future__ import annotations

import csv
import logging
import os
import re
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import duckdb
import numpy as np
import requests
from dotenv import load_dotenv
from flask import Flask, jsonify, make_response, request, send_from_directory
from flask_cors import CORS
from sgp4.api import Satrec, jday

# ── 向量化 SGP4 ───────────────────────────────────────────────────────────────
try:
    from sgp4.api import SatrecArray as _SatrecArray
    _HAS_SATREC_ARRAY = True
except ImportError:
    _HAS_SATREC_ARRAY = False

# ── KD-tree（for 接近事件掃描）────────────────────────────────────────────────
try:
    from scipy.spatial import cKDTree as _cKDTree
    _HAS_KDTREE = True
except ImportError:
    _HAS_KDTREE = False

# ─────────────────────────────────────────────────────────────────────────────
_SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DB  = str(_SCRIPT_DIR / "space_db_slim.duckdb")
RAW_TABLE   = "raw_tle_archive"
META_TABLE  = "sat_n2yo_metadata"

load_dotenv(_SCRIPT_DIR / ".env")

logging.basicConfig(
    level=os.getenv("LOG_LEVEL", "INFO").upper(),
    format="%(asctime)s %(levelname)s %(name)s %(message)s",
)
logger = logging.getLogger("scenario04_adv01")

if _HAS_SATREC_ARRAY:
    logger.info("SatrecArray（向量化 SGP4）已啟用")
else:
    logger.warning("SatrecArray 未可用（sgp4 < 2.0？），退回逐顆傳播")

R_EARTH_KM = 6378.137
F_EARTH    = 1 / 298.257223563
E2         = F_EARTH * (2 - F_EARTH)

DB_PATH               = Path(os.getenv("DB_PATH", DEFAULT_DB))
SAT_META_CSV          = _SCRIPT_DIR / "sat_metadata.csv"
HOST                  = os.getenv("HOST", "0.0.0.0")
PORT                  = int(os.getenv("PORT", "5011"))
STATS_TTL             = int(os.getenv("STATS_TTL", "600"))
_CONJ_THRESHOLD_KM    = float(os.getenv("CONJ_THRESHOLD_KM", "10.0"))
_CONJ_TTL             = int(os.getenv("CONJ_TTL", "120"))

# ── 離線資源 ──────────────────────────────────────────────────────────────────
_GLOBE_TEXTURE_LOCAL = _SCRIPT_DIR / "data" / "globe_texture.jpg"
_CESIUM_LOCAL_DIR    = _SCRIPT_DIR / "data" / "cesium"
_BORDERS_LOCAL       = _SCRIPT_DIR / "data" / "borders.geojson"
_NE_BORDERS_URL      = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector"
    "/master/geojson/ne_110m_admin_0_countries.geojson"
)
_globe_texture_cache: bytes | None = None
_borders_cache:       bytes | None = None


def _get_globe_texture() -> bytes | None:
    global _globe_texture_cache
    if _globe_texture_cache is not None:
        return _globe_texture_cache
    if _GLOBE_TEXTURE_LOCAL.exists():
        try:
            _globe_texture_cache = _GLOBE_TEXTURE_LOCAL.read_bytes()
            logger.info("Globe 貼圖載入: %d bytes", len(_globe_texture_cache))
            return _globe_texture_cache
        except Exception as exc:
            logger.warning("Globe 貼圖讀取失敗: %s", exc)
    logger.error("Globe 貼圖不存在: %s", _GLOBE_TEXTURE_LOCAL)
    return None


# ─────────────────────────────────────────────────────────────────────────────
# 分類函式
# ─────────────────────────────────────────────────────────────────────────────
_COUNTRY_MAP: dict[str, str] = {
    "United States":                  "美國",
    "Commonwealth of Independent States": "俄羅斯/蘇聯",
    "People's Republic of China":     "中國",
    "United Kingdom":                 "英國",
    "France":                         "法國",
    "Japan":                          "日本",
    "India":                          "印度",
    "European Space Agency":          "ESA",
    "Italy":                          "義大利",
    "Canada":                         "加拿大",
    "Germany":                        "德國",
    "South Korea":                    "韓國",
    "Israel":                         "以色列",
    "Australia":                      "澳洲",
    "Brazil":                         "巴西",
    "Argentina":                      "阿根廷",
    "United Arab Emirates":           "阿聯",
    "Iran":                           "伊朗",
    "North Korea":                    "朝鮮",
    "Globalstar":                     "Globalstar",
    "ORBCOMM":                        "ORBCOMM",
    "SOCIETE EUROPEENNE":             "SES",
    "International Telecommunications": "ITSO/Intelsat",
    "EUTELSAT":                       "Eutelsat",
    "TBD":                            "不明",
}
_KNOWN_LABELS = set(_COUNTRY_MAP.values())


def classify_country(source_code: str | None) -> str:
    if not source_code:
        return "不明"
    sc = source_code.strip()
    if sc in _KNOWN_LABELS:
        return sc
    for key, label in _COUNTRY_MAP.items():
        if key.lower() in sc.lower():
            return label
    return "其他"


_CONSTELLATION_RULES: list[tuple[str, list[str]]] = [
    ("Starlink",     ["STARLINK"]),
    ("OneWeb",       ["ONEWEB"]),
    ("Kuiper",       ["KUIPER"]),
    ("千帆/Qianfan", ["QIANFAN", "SPACESAIL"]),
    ("Iridium",      ["IRIDIUM"]),
    ("Globalstar",   ["GLOBALSTAR"]),
    ("Orbcomm",      ["ORBCOMM"]),
    ("Planet/Flock", ["FLOCK", "DOVE", "SKYSAT"]),
    ("Spire",        ["LEMUR", "SPIRE"]),
    ("Telesat LEO",  ["TELESAT"]),
    ("吉林/Jilin",   ["JILIN"]),
    ("遙感/Yaogan",  ["YAOGAN"]),
    ("高分",         ["GAOFEN"]),
    ("風雲",         ["FENGYUN", "FY-"]),
]


def classify_constellation(name: str) -> str:
    n = name.upper()
    for label, kws in _CONSTELLATION_RULES:
        if any(k in n for k in kws):
            return label
    return None  # type: ignore[return-value]


def classify_purpose(name: str) -> str:
    n = name.upper()
    if " DEB" in n or n.endswith(" DEB") or "DEBRIS" in n:
        return "碎片"
    if " R/B" in n or n.endswith(" R/B") or " RB" in n or "ROCKET BODY" in n:
        return "火箭體"
    if "OBJECT" in n:
        return "不明物體"
    return "有效載荷"


def classify_era(launch_date: datetime | None, intl_code: str | None) -> str:
    if launch_date is None and intl_code:
        m = re.match(r"^(\d{4})", str(intl_code))
        if m:
            try:
                launch_date = datetime(int(m.group(1)), 1, 1, tzinfo=timezone.utc)
            except ValueError:
                pass
    if launch_date is None:
        return "不明"
    if launch_date.tzinfo is None:
        launch_date = launch_date.replace(tzinfo=timezone.utc)
    delta_days = (datetime.now(timezone.utc) - launch_date).days
    if delta_days < 365:
        return "< 1 年"
    if delta_days < 365 * 5:
        return "1–5 年"
    if delta_days < 365 * 10:
        return "5–10 年"
    return "> 10 年"


# ─────────────────────────────────────────────────────────────────────────────
# 衛星索引與統計快取
# ─────────────────────────────────────────────────────────────────────────────
_stats_cache:     dict[str, Any] = {}
_stats_loaded_at: float = 0.0
_sat_index:       dict[int, dict[str, Any]] = {}
_index_loaded_at: float = 0.0
_INDEX_TTL = 600


def load_sat_metadata_csv() -> dict[int, dict[str, str]]:
    if not SAT_META_CSV.exists():
        return {}
    result: dict[int, dict[str, str]] = {}
    try:
        with SAT_META_CSV.open(encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f):
                raw_id = row.get("norad_id", "").strip()
                if not raw_id:
                    continue
                try:
                    nid = int(raw_id)
                except ValueError:
                    continue
                result[nid] = {k: (v.strip() if v else "") for k, v in row.items()
                               if k != "norad_id"}
        logger.info("sat_metadata.csv: %d 筆", len(result))
    except Exception as exc:
        logger.error("sat_metadata.csv 讀取失敗: %s", exc)
    return result


def _resolve_db() -> Path | None:
    if DB_PATH.exists():
        return DB_PATH
    # 若指定路徑不存在，嘗試同目錄下的備用 DB
    alt_name = "space_db.duckdb" if "slim" in DB_PATH.name else "space_db_slim.duckdb"
    alt = DB_PATH.parent / alt_name
    if alt.exists():
        logger.warning("DB %s 不存在，改用 %s", DB_PATH.name, alt.name)
        return alt
    logger.error("找不到資料庫: %s", DB_PATH)
    return None


def build_sat_index() -> dict[int, dict[str, Any]]:
    db = _resolve_db()
    if db is None:
        return {}
    logger.info("建立衛星索引中…")
    t0 = time.monotonic()
    try:
        with duckdb.connect(str(db), read_only=True) as con:
            rows = con.execute(f"""
                SELECT
                    r.norad_id, r.object_name, r.line1, r.line2,
                    m.source_code, m.launch_date, m.intl_code
                FROM {RAW_TABLE} r
                LEFT JOIN {META_TABLE} m ON r.norad_id = m.norad_id
                WHERE r.line1 IS NOT NULL AND r.line2 IS NOT NULL
                QUALIFY ROW_NUMBER() OVER (
                    PARTITION BY r.norad_id ORDER BY r.epoch_utc DESC
                ) = 1
            """).fetchall()
    except Exception as exc:
        logger.error("建立索引失敗: %s", exc)
        return {}

    csv_meta = load_sat_metadata_csv()
    idx: dict[int, dict[str, Any]] = {}
    for norad_id, raw_name, l1, l2, db_src, db_launch, db_intl in rows:
        nid  = int(norad_id)
        name = (raw_name or "").strip().lstrip("0 ") or f"OBJECT {nid}"
        ov   = csv_meta.get(nid, {})
        final_name = ov.get("name_en")    or name
        final_src  = ov.get("source_code") or db_src
        final_intl = ov.get("intl_code")  or db_intl
        csv_date_str = ov.get("launch_date", "")
        if csv_date_str:
            try:
                final_launch: datetime | None = datetime.strptime(
                    csv_date_str, "%Y-%m-%d").replace(tzinfo=timezone.utc)
            except ValueError:
                final_launch = db_launch
        else:
            final_launch = db_launch
        csv_purpose = ov.get("purpose", "")
        csv_constel = ov.get("constellation", "")
        idx[nid] = {
            "name":          final_name,
            "line1":         l1.strip() if l1 else "",
            "line2":         l2.strip() if l2 else "",
            "country":       classify_country(final_src),
            "purpose":       csv_purpose if csv_purpose else classify_purpose(final_name),
            "era":           classify_era(final_launch, final_intl),
            "constellation": csv_constel if csv_constel else classify_constellation(final_name),
        }
    elapsed = time.monotonic() - t0
    logger.info("衛星索引完成: %d 筆，耗時 %.1f s", len(idx), elapsed)
    return idx


def get_sat_index() -> dict[int, dict[str, Any]]:
    global _sat_index, _index_loaded_at
    if not _sat_index or (time.monotonic() - _index_loaded_at) > _INDEX_TTL:
        _sat_index = build_sat_index()
        _index_loaded_at = time.monotonic()
    return _sat_index


def build_stats(idx: dict[int, dict[str, Any]]) -> dict[str, Any]:
    country: dict[str, int] = {}
    purpose: dict[str, int] = {}
    era:     dict[str, int] = {}
    constel: dict[str, int] = {}
    for info in idx.values():
        country[info["country"]] = country.get(info["country"], 0) + 1
        purpose[info["purpose"]] = purpose.get(info["purpose"], 0) + 1
        era[info["era"]]         = era.get(info["era"], 0) + 1
        c = info["constellation"] or "其他衛星"
        constel[c] = constel.get(c, 0) + 1

    def _sorted(d: dict[str, int]) -> list[dict]:
        return [{"label": k, "count": v} for k, v in sorted(d.items(), key=lambda x: -x[1])]

    era_order = ["< 1 年", "1–5 年", "5–10 年", "> 10 年", "不明"]
    era_sorted = sorted(era.items(),
                        key=lambda x: era_order.index(x[0]) if x[0] in era_order else 99)
    return {
        "total":         len(idx),
        "country":       _sorted(country),
        "purpose":       _sorted(purpose),
        "era":           [{"label": k, "count": v} for k, v in era_sorted],
        "constellation": _sorted(constel),
        "updated_at":    datetime.now(timezone.utc).isoformat(),
    }


def get_stats() -> dict[str, Any]:
    global _stats_cache, _stats_loaded_at
    if not _stats_cache or (time.monotonic() - _stats_loaded_at) > STATS_TTL:
        idx = get_sat_index()
        _stats_cache = build_stats(idx)
        _stats_loaded_at = time.monotonic()
    return _stats_cache


# ─────────────────────────────────────────────────────────────────────────────
# 1-A  向量化 SGP4
# ─────────────────────────────────────────────────────────────────────────────

def eci_to_llh_batch(r_arr: np.ndarray, t: datetime) -> np.ndarray:
    """
    向量化 ECI→geodetic 轉換。
    r_arr : (N, 3) ECI 位置 km
    returns: (N, 3)  [lat_deg, lon_deg, alt_km]
    """
    x, y, z = r_arr[:, 0], r_arr[:, 1], r_arr[:, 2]
    jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute,
                  t.second + t.microsecond * 1e-6)
    T_cent  = ((jd - 2451545.0) + fr) / 36525.0
    gmst    = np.deg2rad(
        (280.46061837 + 360.98564736629 * (jd - 2451545.0 + fr)
         + 0.000387933 * T_cent ** 2) % 360.0)
    xe = np.cos(gmst) * x + np.sin(gmst) * y
    ye = -np.sin(gmst) * x + np.cos(gmst) * y
    ze = z
    lon = np.arctan2(ye, xe)
    rr  = np.sqrt(xe ** 2 + ye ** 2)
    lat = np.arctan2(ze, rr * (1.0 - E2))
    alt = np.zeros(len(r_arr), dtype=float)
    for _ in range(5):
        sl    = np.sin(lat)
        N_arr = R_EARTH_KM / np.sqrt(1.0 - E2 * sl ** 2)
        cl    = np.cos(lat)
        alt   = np.where(np.abs(cl) > 1e-9,
                         rr / cl - N_arr,
                         np.abs(ze) / (1.0 - E2) - N_arr)
        lat   = np.arctan2(ze, rr * (1.0 - E2 * (N_arr / (N_arr + alt))))
    return np.column_stack([np.rad2deg(lat), np.rad2deg(lon), alt])


def _sgp4_propagate_raw(
    nids:   list[int],
    line1s: list[str],
    line2s: list[str],
    t:      datetime,
) -> tuple[np.ndarray, np.ndarray]:
    """
    內部：批次 SGP4 傳播。
    Returns
      err_arr : (N,) int  — 0 = success
      r_arr   : (N, 3) float  — ECI km（err != 0 的列為 0）
    """
    n   = len(nids)
    jd0, fr0 = jday(t.year, t.month, t.day, t.hour, t.minute,
                    t.second + t.microsecond * 1e-6)

    if _HAS_SATREC_ARRAY and n > 1:
        try:
            sats = _SatrecArray.twoline2rv(line1s, line2s)
            e_raw, r_raw, _ = sats.sgp4(np.array([jd0]), np.array([fr0]))
            # e_raw: (N,1)  r_raw: (N,1,3)
            return e_raw[:, 0].astype(int), r_raw[:, 0, :]
        except Exception as exc:
            logger.debug("SatrecArray 傳播失敗，退回逐顆: %s", exc)

    # Sequential fallback
    err_arr = np.zeros(n, dtype=int)
    r_arr   = np.zeros((n, 3), dtype=float)
    for i, (l1, l2) in enumerate(zip(line1s, line2s)):
        try:
            sat = Satrec.twoline2rv(l1, l2)
            err, r, _ = sat.sgp4(jd0, fr0)
            err_arr[i] = err
            if err == 0:
                r_arr[i] = r
        except Exception:
            err_arr[i] = 1
    return err_arr, r_arr


def propagate_batch(
    nids: list[int],
    idx:  dict[int, dict[str, Any]],
) -> list[tuple[float, float, float] | None]:
    """
    公開 API：批次傳播，回傳 (lat, lon, alt_km) 或 None（每顆對應）。
    """
    if not nids:
        return []
    line1s = [idx[n]["line1"] for n in nids]
    line2s = [idx[n]["line2"] for n in nids]
    t = datetime.now(timezone.utc)

    err_arr, r_arr = _sgp4_propagate_raw(nids, line1s, line2s, t)
    valid = err_arr == 0

    if not valid.any():
        return [None] * len(nids)

    llh = np.full((len(nids), 3), np.nan)
    llh[valid] = eci_to_llh_batch(r_arr[valid], t)

    results: list[tuple[float, float, float] | None] = []
    for i in range(len(nids)):
        if not valid[i] or np.isnan(llh[i, 0]):
            results.append(None)
            continue
        lat, lon, alt = float(llh[i, 0]), float(llh[i, 1]), float(llh[i, 2])
        if not (-500.0 < alt < 80_000.0):
            results.append(None)
            continue
        results.append((lat, lon, alt))
    return results


def propagate_now(line1: str, line2: str) -> tuple[float, float, float] | None:
    """單顆傳播（保留，供搜尋/單查路由使用）。"""
    try:
        sat = Satrec.twoline2rv(line1, line2)
        t   = datetime.now(timezone.utc)
        jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute,
                      t.second + t.microsecond * 1e-6)
        err, r_eci, _ = sat.sgp4(jd, fr)
        if err != 0:
            return None
        llh = eci_to_llh_batch(np.array([r_eci], dtype=float), t)
        lat, lon, alt = float(llh[0, 0]), float(llh[0, 1]), float(llh[0, 2])
        if not (-500.0 < alt < 80_000.0):
            return None
        return lat, lon, alt
    except Exception:
        return None


# ─────────────────────────────────────────────────────────────────────────────
# 2-A  即時近距離配對掃描
# ─────────────────────────────────────────────────────────────────────────────
_conj_cache:     dict[str, Any] | None = None
_conj_loaded_at: float = 0.0


def build_conjunction_summary(
    threshold_km: float = _CONJ_THRESHOLD_KM,
    max_pairs:    int   = 200,
) -> dict[str, Any]:
    """
    對全目錄進行向量化 SGP4 傳播，在 ECI 空間建 cKDTree，
    query_pairs 找出當前位置距離 ≤ threshold_km 的衛星配對。

    注意：這是「當前瞬間位置近距離」，不含 TCA/Pc 計算。
    詳細接近分析請使用 conjunction_pipeline.py / ca_pipeline_kdtree_v2_fixed.py。
    """
    if not _HAS_KDTREE:
        return {"error": "scipy 未安裝，無法執行接近事件掃描",
                "count": 0, "pairs": [], "threshold_km": threshold_km}

    idx = get_sat_index()
    if len(idx) < 2:
        return {"count": 0, "pairs": [], "threshold_km": threshold_km, "total_scanned": 0}

    all_nids = list(idx.keys())
    line1s   = [idx[n]["line1"] for n in all_nids]
    line2s   = [idx[n]["line2"] for n in all_nids]
    t = datetime.now(timezone.utc)

    t0 = time.monotonic()
    err_arr, r_arr = _sgp4_propagate_raw(all_nids, line1s, line2s, t)
    t_sgp4 = time.monotonic() - t0

    # 過濾 SGP4 成功的衛星
    ok = (err_arr == 0)
    ok_nids = [all_nids[i] for i in range(len(all_nids)) if ok[i]]
    ok_r    = r_arr[ok]                        # (M, 3)

    if len(ok_nids) < 2:
        return {"count": 0, "pairs": [], "threshold_km": threshold_km,
                "total_scanned": len(ok_nids)}

    # 批次 ECI→LLH，過濾合理高度
    ok_llh  = eci_to_llh_batch(ok_r, t)       # (M, 3)
    alt_ok  = (ok_llh[:, 2] > -500.0) & (ok_llh[:, 2] < 80_000.0)
    filt_nids = [ok_nids[i] for i in range(len(ok_nids)) if alt_ok[i]]
    filt_r    = ok_r[alt_ok]                   # (F, 3) ECI
    filt_llh  = ok_llh[alt_ok]                # (F, 3) LLH

    if len(filt_nids) < 2:
        return {"count": 0, "pairs": [], "threshold_km": threshold_km,
                "total_scanned": len(filt_nids)}

    t1 = time.monotonic()
    tree      = _cKDTree(filt_r)
    pairs_set = tree.query_pairs(threshold_km)
    t_kd      = time.monotonic() - t1

    logger.info(
        "接近事件掃描: %d 有效衛星，閾值 %.0f km，配對 %d，SGP4 %.2f s，KD %.2f s",
        len(filt_nids), threshold_km, len(pairs_set), t_sgp4, t_kd,
    )

    pairs: list[dict[str, Any]] = []
    for i, j in pairs_set:
        miss = float(np.linalg.norm(filt_r[i] - filt_r[j]))
        nid_a, nid_b = filt_nids[i], filt_nids[j]
        pairs.append({
            "primary_norad":    nid_a,
            "primary_name":     idx[nid_a]["name"],
            "primary_purpose":  idx[nid_a]["purpose"],
            "primary_lat":      round(float(filt_llh[i, 0]), 4),
            "primary_lon":      round(float(filt_llh[i, 1]), 4),
            "primary_alt_km":   round(float(filt_llh[i, 2]), 1),
            "secondary_norad":  nid_b,
            "secondary_name":   idx[nid_b]["name"],
            "secondary_purpose":idx[nid_b]["purpose"],
            "secondary_lat":    round(float(filt_llh[j, 0]), 4),
            "secondary_lon":    round(float(filt_llh[j, 1]), 4),
            "secondary_alt_km": round(float(filt_llh[j, 2]), 1),
            "miss_km":          round(miss, 3),
        })

    pairs.sort(key=lambda x: x["miss_km"])
    if len(pairs) > max_pairs:
        pairs = pairs[:max_pairs]

    elapsed = time.monotonic() - t0
    return {
        "count":         len(pairs),
        "threshold_km":  threshold_km,
        "max_pairs":     max_pairs,
        "pairs":         pairs,
        "total_scanned": len(filt_nids),
        "elapsed_sec":   round(elapsed, 2),
        "vectorized":    _HAS_SATREC_ARRAY,
        "timestamp":     datetime.now(timezone.utc).isoformat(),
    }


def get_conjunctions(threshold_km: float = _CONJ_THRESHOLD_KM) -> dict[str, Any]:
    global _conj_cache, _conj_loaded_at
    if _conj_cache is None or (time.monotonic() - _conj_loaded_at) > _CONJ_TTL:
        _conj_cache = build_conjunction_summary(threshold_km=threshold_km)
        _conj_loaded_at = time.monotonic()
    return _conj_cache


# ─────────────────────────────────────────────────────────────────────────────
# 顏色映射
# ─────────────────────────────────────────────────────────────────────────────
_PURPOSE_COLORS = {
    "有效載荷": "#4CAF50", "碎片": "#FF9800",
    "火箭體":   "#9E9E9E", "不明物體": "#607D8B",
}
_COUNTRY_COLORS = {
    "美國": "#3F51B5", "俄羅斯/蘇聯": "#F44336", "中國": "#FF5722",
    "英國": "#2196F3", "法國": "#9C27B0", "日本": "#E91E63",
    "印度": "#FF9800", "ESA": "#00BCD4", "其他": "#78909C", "不明": "#455A64",
}
_CONSTELLATION_COLORS = {
    "Starlink": "#1565C0", "OneWeb": "#00897B", "Kuiper": "#FF8F00",
    "千帆/Qianfan": "#C62828", "Iridium": "#558B2F", "Globalstar": "#6A1B9A",
    "Planet/Flock": "#2E7D32", "Spire": "#00838F", "吉林/Jilin": "#AD1457",
    "遙感/Yaogan": "#B71C1C", "高分": "#E64A19", "風雲": "#0277BD",
    "其他衛星": "#546E7A",
}
_ERA_COLORS = {
    "< 1 年": "#F44336", "1–5 年": "#FF9800",
    "5–10 年": "#4CAF50", "> 10 年": "#607D8B", "不明": "#455A64",
}


def get_color(ftype: str, label: str) -> str:
    maps = {
        "purpose":       _PURPOSE_COLORS,
        "country":       _COUNTRY_COLORS,
        "constellation": _CONSTELLATION_COLORS,
        "era":           _ERA_COLORS,
    }
    return maps.get(ftype, {}).get(label, "#78909C")


# ─────────────────────────────────────────────────────────────────────────────
# SSN 地面觀測站（嵌入資料）
# ─────────────────────────────────────────────────────────────────────────────
_SSN_STATIONS_GEOJSON = {
    "type": "FeatureCollection",
    "features": [
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-106.6599,33.8172]},"properties":{"name":"GEODSS Socorro","type":"光學/電光","location":"White Sands Missile Range, NM, USA","status":"active","notes":"深空目標光學追蹤"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-156.2578,20.7088]},"properties":{"name":"GEODSS Maui (AMOS)","type":"光學/電光","location":"Haleakalā, Hawaii, USA","status":"active","notes":"深空目標光學追蹤"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[72.4522,-7.4117]},"properties":{"name":"GEODSS Diego Garcia","type":"光學/電光","location":"British Indian Ocean Territory","status":"active","notes":"深空目標光學追蹤"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-156.2600,20.7100]},"properties":{"name":"MAUI Space Surveillance (MSSS)","type":"光學/電光","location":"Haleakalā, Hawaii, USA","status":"active","notes":"先進光電系統，多光譜成像"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-156.2565,20.7076]},"properties":{"name":"AEOS Telescope","type":"光學/電光","location":"Haleakalā, Hawaii, USA","status":"active","notes":"先進電光感測器（3.67m 望遠鏡）"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-14.40,-7.97]},"properties":{"name":"Ascension Range Radar","type":"雷達","location":"Ascension Island","status":"active","notes":"南大西洋遙測/追蹤站"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-86.21,30.57]},"properties":{"name":"AN/FPS-85","type":"雷達","location":"Eglin AFB, Florida, USA","status":"active","notes":"SSN最大功率相控陣雷達；第20太空監視中隊操作"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[167.750,8.730]},"properties":{"name":"Space Fence (AN/FPS-133)","type":"雷達","location":"Kwajalein Atoll, Marshall Islands","status":"active","notes":"S波段相控陣雷達，2020年起作戰；可追蹤10cm以下碎片"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[174.14,52.72]},"properties":{"name":"AN/FPS-108 Cobra Dane","type":"雷達","location":"Shemya Island, Alaska, USA","status":"active","notes":"相控陣雷達，兼飛彈預警與太空追蹤"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[31.13,70.37]},"properties":{"name":"GLOBUS II","type":"雷達","location":"Vardø, Norway","status":"active","notes":"X波段碟形雷達；挪威情報局（NIS）操作，數據納入SSN"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-71.518,42.261]},"properties":{"name":"HUSIR (Haystack)","type":"雷達","location":"Westford, Massachusetts, USA","status":"active","notes":"超寬頻雷達，MIT林肯實驗室"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-71.87,42.28]},"properties":{"name":"Millstone Hill Radar","type":"雷達","location":"North Grafton, Massachusetts, USA","status":"active","notes":"MIT林肯實驗室追蹤雷達"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[167.755,8.735]},"properties":{"name":"ALTAIR","type":"雷達","location":"Kwajalein Atoll, Marshall Islands","status":"active","notes":"Reagan Test Site (RTS) 深空雷達"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-96.37,47.12]},"properties":{"name":"PARCS (AN/FPQ-16)","type":"雷達","location":"Cavalier AFS, North Dakota, USA","status":"active","notes":"飛彈預警/太空監視相控陣雷達；第10太空預警中隊"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[167.745,8.725]},"properties":{"name":"ALCOR","type":"雷達","location":"Kwajalein Atoll, Marshall Islands","status":"active","notes":"C波段成像雷達，Reagan Test Site (RTS)"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-70.54,41.75]},"properties":{"name":"Cape Cod SFS (PAVE PAWS)","type":"飛彈預警/協作","location":"Bourne, Massachusetts, USA","status":"active","notes":"AN/FPS-123，UHF相控陣；第6太空預警中隊"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-121.35,39.14]},"properties":{"name":"Beale AFB (PAVE PAWS)","type":"飛彈預警/協作","location":"California, USA","status":"active","notes":"AN/FPS-123，UHF相控陣；第7太空預警中隊"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-149.19,64.30]},"properties":{"name":"Clear SFS (BMEWS)","type":"飛彈預警/協作","location":"Alaska, USA","status":"active","notes":"AN/FPS-120，UHF相控陣；第13太空預警中隊"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-0.67,54.36]},"properties":{"name":"RAF Fylingdales (BMEWS)","type":"飛彈預警/協作","location":"England, UK","status":"active","notes":"AN/FPS-126，三面相控陣，360°覆蓋；英國皇家空軍操作"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-68.78,76.53]},"properties":{"name":"Thule/Pituffik SB (BMEWS)","type":"飛彈預警/協作","location":"Greenland","status":"active","notes":"AN/FPS-120；兼為SCN衛星追控站"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-71.43,42.97]},"properties":{"name":"New Boston AFS","type":"衛星追控","location":"New Hampshire, USA","status":"active","notes":"23rd Space Operations Squadron SCN站點"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-120.57,34.75]},"properties":{"name":"Vandenberg SFB (SCN)","type":"衛星追控","location":"California, USA","status":"active","notes":"23rd Space Operations Squadron SCN站點"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-158.13,21.36]},"properties":{"name":"Kaena Point","type":"衛星追控","location":"Oahu, Hawaii, USA","status":"active","notes":"23rd Space Operations Squadron SCN站點"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[144.87,13.58]},"properties":{"name":"Guam Remote Tracking Station","type":"衛星追控","location":"Guam","status":"active","notes":"23rd Space Operations Squadron SCN站點"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-1.08,51.07]},"properties":{"name":"RAF Oakhanger","type":"衛星追控","location":"Hampshire, UK","status":"active","notes":"23rd Space Operations Squadron SCN站點（英國）"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[72.4508,-7.413]},"properties":{"name":"Diego Garcia (SCN)","type":"衛星追控","location":"British Indian Ocean Territory","status":"active","notes":"23rd Space Operations Squadron SCN站點（兼GEODSS）"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-104.80,38.70]},"properties":{"name":"Space Surveillance Center (SSC)","type":"數據中心","location":"Cheyenne Mountain Complex, Colorado, USA","status":"active","notes":"SPACECOM太空監視作戰中心"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[40.04,37.95]},"properties":{"name":"AN/FPS-79 Pirinclik","type":"已除役","location":"Diyarbakır, Turkey","status":"decommissioned","notes":"原SSN輔助雷達，已關閉"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-5.609,37.17]},"properties":{"name":"MOSS (Morón)","type":"已除役","location":"Morón Air Base, Spain","status":"decommissioned","notes":"光學站，1997–2012年運作"}},
    ],
}


# ─────────────────────────────────────────────────────────────────────────────
# 嵌入式前端
# ─────────────────────────────────────────────────────────────────────────────
_HTML = r"""<!DOCTYPE html>
<html lang="zh-TW">
<head>
<meta charset="utf-8"/>
<title>太空態勢儀表板 — 進階版 04-adv01（向量化 SGP4 + 近距離掃描 + 搜尋）</title>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<style>
*{box-sizing:border-box;margin:0;padding:0}
html,body{height:100%;overflow:hidden}
body{font-family:Tahoma,sans-serif;background:#0a0e17;color:#c9d1d9;display:flex;flex-direction:column}

/* ── topbar ── */
#topbar{display:flex;align-items:center;gap:12px;padding:5px 14px;
  background:linear-gradient(90deg,#0d1117 0%,#161b22 100%);
  border-bottom:1px solid #1e3a5f;flex-shrink:0;height:46px;overflow-x:auto}
#topbar h1{font-size:13px;color:#58a6ff;white-space:nowrap;letter-spacing:.5px}
.tcards{display:flex;gap:6px;flex:1}
.tcard{background:#1a2744;border:1px solid #1e3a5f;border-radius:4px;
  padding:3px 10px;text-align:center;min-width:80px;white-space:nowrap}
.tcard.clickable{cursor:pointer;transition:.15s}
.tcard.clickable:hover{border-color:#58a6ff;background:#1e3a5f}
.tcard.active-card{border-color:#58a6ff;background:#1e3a6f}
.tcard .val{font-size:16px;font-weight:bold;color:#f0f6ff}
.tcard .lbl{font-size:9px;color:#8b949e;text-transform:uppercase;letter-spacing:.6px}
.tcard.payload .val{color:#4CAF50}
.tcard.debris  .val{color:#FF9800}
.tcard.rocket  .val{color:#9E9E9E}
.tcard.conj    .val{color:#FF6B6B}
#ts{font-size:10px;color:#484f58;white-space:nowrap;flex-shrink:0}

/* ── main layout ── */
#main{display:flex;flex:1;overflow:hidden}
#panel{width:300px;min-width:260px;background:#0d1117;
  border-right:1px solid #21262d;display:flex;flex-direction:column;overflow:hidden}

/* ── search box ── */
#search-box{display:flex;align-items:center;gap:4px;
  padding:6px 10px 4px;border-bottom:1px solid #21262d;flex-shrink:0}
#search-input{flex:1;background:#161b22;border:1px solid #21262d;border-radius:4px;
  color:#c9d1d9;font-size:11px;padding:5px 8px;outline:none;transition:.15s}
#search-input:focus{border-color:#58a6ff}
#search-input::placeholder{color:#484f58}
#search-clear{background:none;border:none;color:#484f58;cursor:pointer;
  font-size:15px;line-height:1;padding:0 2px;transition:.12s;display:none}
#search-clear:hover{color:#c9d1d9}

/* ── tabs ── */
#tabs{display:flex;border-bottom:1px solid #21262d;flex-shrink:0}
.tab{flex:1;padding:6px 0;font-size:11px;color:#8b949e;text-align:center;
  cursor:pointer;border:none;background:transparent;transition:.15s;border-bottom:2px solid transparent}
.tab:hover{color:#c9d1d9;background:#161b22}
.tab.active{color:#58a6ff;border-bottom:2px solid #58a6ff}
#panel-back{display:none;padding:5px 10px 4px;border-bottom:1px solid #21262d;flex-shrink:0}
#panel-back button{background:none;border:none;color:#8b949e;font-size:11px;cursor:pointer;
  display:flex;align-items:center;gap:5px}
#panel-back button:hover{color:#c9d1d9}

/* ── panel body ── */
#panel-body{flex:1;overflow-y:auto;padding:4px 0}
#panel-body::-webkit-scrollbar{width:3px}
#panel-body::-webkit-scrollbar-thumb{background:#21262d;border-radius:2px}

.stat-row{display:flex;align-items:center;gap:7px;padding:6px 12px;
  cursor:pointer;border:none;background:transparent;width:100%;text-align:left;
  border-left:3px solid transparent;transition:.12s}
.stat-row:hover{background:#1c2128}
.stat-row.active{background:#1a2744;border-left-color:#58a6ff}
.stat-dot{width:8px;height:8px;border-radius:50%;flex-shrink:0}
.stat-label{flex:1;font-size:12px;color:#c9d1d9;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.stat-cnt{font-size:11px;color:#58a6ff;font-weight:bold;min-width:40px;text-align:right}
.stat-bar-wrap{width:100%;height:2px;background:#21262d;border-radius:1px}
.stat-bar{height:2px;border-radius:1px;transition:width .3s}

/* ── conjunction rows ── */
.conj-row{padding:6px 10px;border-bottom:1px solid #161b22;cursor:pointer;transition:.12s}
.conj-row:hover{background:#1c2128}
.conj-km{font-size:11px;font-weight:bold;color:#FF6B6B;min-width:50px}
.conj-names{flex:1;font-size:11px;color:#c9d1d9;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.conj-note{font-size:9px;color:#484f58;padding:4px 10px 8px;font-style:italic}

/* ── filter status & layer bar ── */
#filter-status{padding:5px 10px;font-size:11px;color:#8b949e;
  border-top:1px solid #21262d;flex-shrink:0;min-height:26px}
#layer-bar{border-top:1px solid #21262d;background:#0d1117;padding:5px 0 3px;flex-shrink:0}
#layer-bar .sec{padding:3px 12px 2px;font-size:10px;color:#484f58;text-transform:uppercase;letter-spacing:.8px}
#layer-bar label{display:flex;align-items:center;gap:7px;padding:4px 12px;cursor:pointer;
  color:#8b949e;font-size:11px;transition:.12s;user-select:none}
#layer-bar label:hover{background:#1c2128;color:#c9d1d9}
#layer-bar input[type=checkbox]{margin:0;cursor:pointer;accent-color:#58a6ff}

/* ── cesium ── */
#cesiumContainer{flex:1;position:relative}
#warn{position:absolute;top:10px;left:50%;transform:translateX(-50%);
  color:#fff;padding:5px 14px;border-radius:4px;font-size:11px;z-index:10;
  pointer-events:none;display:none;white-space:nowrap}
#loading{position:absolute;inset:0;background:rgba(10,14,23,.75);
  display:flex;align-items:center;justify-content:center;z-index:20;
  font-size:13px;color:#58a6ff;letter-spacing:1px}
.cesium-widget-credits{display:none!important}
.cesium-infoBox-description{background:#1e1e1e!important;color:#e8e8e8!important}
.cesium-infoBox-description a{color:#8ab4f8!important}
</style>
<!-- 本機 Cesium 載入（data/cesium/ = CesiumJS 1.114） -->
<script>
(function(){
  function _css(h){var l=document.createElement('link');l.rel='stylesheet';l.href=h;document.head.appendChild(l);}
  function _js(src,ok,fail){var s=document.createElement('script');s.src=src;s.onload=ok;s.onerror=fail||function(){};document.head.appendChild(s);}
  function _err(m){
    function _s(){var el=document.getElementById('warn');if(el){el.textContent=m;el.style.background='rgba(244,67,54,.88)';el.style.display='block';}}
    if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',_s);else _s();
  }
  _css('/cesium/Widgets/widgets.css');
  _js('/cesium/Cesium.js',function(){startApp();},function(){_err('本機 Cesium 未安裝，請將 CesiumJS 解壓至 data/cesium/');});
}());
</script>
</head>
<body>

<div id="topbar">
  <h1>&#127760; 太空態勢儀表板（進階版）</h1>
  <div class="tcards">
    <div class="tcard"><div class="val" id="c-total">—</div><div class="lbl">追蹤物件</div></div>
    <div class="tcard payload"><div class="val" id="c-payload">—</div><div class="lbl">有效載荷</div></div>
    <div class="tcard debris"><div class="val" id="c-debris">—</div><div class="lbl">碎片</div></div>
    <div class="tcard rocket"><div class="val" id="c-rocket">—</div><div class="lbl">火箭體</div></div>
    <div class="tcard"><div class="val" id="c-new">—</div><div class="lbl">近 1 年</div></div>
    <div class="tcard conj clickable" id="conj-card" onclick="toggleConjPanel()" title="點選查看近距離配對清單">
      <div class="val" id="c-conj">...</div>
      <div class="lbl">近距離配對</div>
    </div>
  </div>
  <div id="ts">—</div>
</div>

<div id="main">
  <div id="panel">
    <!-- 4-A 搜尋框 -->
    <div id="search-box">
      <input id="search-input" type="text" placeholder="搜尋 NORAD ID 或衛星名稱..." autocomplete="off"/>
      <button id="search-clear" onclick="clearSearch()">×</button>
    </div>
    <!-- 面板標題列（分頁 or 返回按鈕） -->
    <div id="tabs">
      <button class="tab active" data-tab="country"       onclick="switchTab('country')">國家</button>
      <button class="tab"        data-tab="purpose"       onclick="switchTab('purpose')">用途</button>
      <button class="tab"        data-tab="era"           onclick="switchTab('era')">年代</button>
      <button class="tab"        data-tab="constellation" onclick="switchTab('constellation')">星座</button>
    </div>
    <div id="panel-back">
      <button onclick="backToTabs()">&#8592; 返回</button>
    </div>
    <div id="panel-body"></div>
    <div id="filter-status">點選左側分類以在地圖上顯示（全量，無抽樣）</div>
    <div id="layer-bar">
      <div class="sec">向量圖層</div>
      <label>
        <input type="checkbox" id="chk-borders" checked onchange="toggleLayer('borders',this)"/>
        <span style="width:9px;height:9px;border-radius:50%;background:#FFD600;display:inline-block;flex-shrink:0"></span>
        <span>全球國界</span>
      </label>
      <label>
        <input type="checkbox" id="chk-ssn" onchange="toggleLayer('ssn',this)"/>
        <span style="width:9px;height:9px;border-radius:50%;background:#00E5FF;display:inline-block;flex-shrink:0"></span>
        <span>SSN 地面觀測站</span>
      </label>
    </div>
  </div>

  <div id="cesiumContainer">
    <div id="warn"></div>
    <div id="loading">正在初始化地球...</div>
  </div>
</div>

<script>
'use strict';
function startApp(){

let viewer=null, satDs=null;
let statsData=null, conjData=null;
let activeTab='country', activeFtype=null, activeFval=null;
// panelMode: 'tabs' | 'search' | 'conjunctions'
let panelMode='tabs';
const entMap=new Map();
let borderDs=null, ssnDs=null;
let _searchTimer=null;

const SSN_TYPE_COLORS={
  '光學/電光':'#FFC107','雷達':'#00E5FF','飛彈預警/協作':'#FF7043',
  '衛星追控':'#66BB6A','數據中心':'#CE93D8','已除役':'#78909C',
};
const PURPOSE_C  ={有效載荷:'#4CAF50',碎片:'#FF9800',火箭體:'#9E9E9E',不明物體:'#607D8B'};
const COUNTRY_C  ={美國:'#3F51B5','俄羅斯/蘇聯':'#F44336',中國:'#FF5722',英國:'#2196F3',
                   法國:'#9C27B0',日本:'#E91E63',印度:'#FF9800',ESA:'#00BCD4',
                   其他:'#78909C',不明:'#455A64'};
const CONSTEL_C  ={Starlink:'#1565C0',OneWeb:'#00897B',Kuiper:'#FF8F00',
                   '千帆/Qianfan':'#C62828',Iridium:'#558B2F',Globalstar:'#6A1B9A',
                   'Planet/Flock':'#2E7D32',Spire:'#00838F','吉林/Jilin':'#AD1457',
                   '遙感/Yaogan':'#B71C1C',高分:'#E64A19',風雲:'#0277BD',其他衛星:'#546E7A'};
const ERA_C      ={'< 1 年':'#F44336','1–5 年':'#FF9800','5–10 年':'#4CAF50',
                   '> 10 年':'#607D8B',不明:'#455A64'};

function getColor(ftype,label){
  const m={country:COUNTRY_C,purpose:PURPOSE_C,era:ERA_C,constellation:CONSTEL_C};
  return (m[ftype]||{})[label]||'#78909C';
}

// ── Cesium 初始化 ──────────────────────────────────────────────────────────
async function initCesium(){
  Cesium.Ion.defaultAccessToken='';
  try{ Cesium.IonResource.fromAssetId=function(){return Promise.reject(new Error('Ion disabled'));}; }catch(e){}
  const opts={
    animation:false, timeline:false, baseLayerPicker:false,
    imageryProvider:new Cesium.SingleTileImageryProvider({
      url:'/api/globe_texture',
      rectangle:Cesium.Rectangle.fromDegrees(-180,-90,180,90),
      credit:'NASA Blue Marble © NASA Earth Observatory',
    }),
    terrainProvider:new Cesium.EllipsoidTerrainProvider(),
    sceneModePicker:true, infoBox:true, geocoder:false,
    homeButton:true, navigationHelpButton:false, selectionIndicator:true,
  };
  viewer=new Cesium.Viewer('cesiumContainer',opts);
  viewer.cesiumWidget.creditContainer.style.display='none';
  viewer.scene.globe.enableLighting=true;
  viewer.scene.globe.depthTestAgainstTerrain=false;
  satDs=new Cesium.CustomDataSource('satellites');
  await viewer.dataSources.add(satDs);
  document.getElementById('loading').style.display='none';
}

// ── 全球國界 ────────────────────────────────────────────────────────────────
async function loadBordersLayer(){
  try{
    const ds=await Cesium.GeoJsonDataSource.load('/api/layers/borders');
    const lineColor=Cesium.Color.fromCssColorString('#FFD600').withAlpha(0.85);
    [...ds.entities.values].forEach(ent=>{
      if(ent.polygon){
        const hier=ent.polygon.hierarchy&&ent.polygon.hierarchy.getValue(Cesium.JulianDate.now());
        if(hier&&hier.positions&&hier.positions.length){
          ds.entities.add({polyline:{
            positions:[...hier.positions,hier.positions[0]],
            width:1.5, clampToGround:true,
            material:new Cesium.ColorMaterialProperty(lineColor),
            arcType:Cesium.ArcType.GEODESIC,
          }});
        }
        ent.polygon.show=new Cesium.ConstantProperty(false);
      }
      if(ent.label) ent.label.show=new Cesium.ConstantProperty(false);
    });
    borderDs=ds;
    await viewer.dataSources.add(ds);
  }catch(e){
    console.warn('全球國界載入失敗',e);
    const chk=document.getElementById('chk-borders');
    if(chk) chk.checked=false;
  }
}

// ── 統計 ──────────────────────────────────────────────────────────────────
async function loadStats(){
  try{
    const r=await fetch('/api/stats');
    if(!r.ok) throw new Error('HTTP '+r.status);
    statsData=await r.json();
    renderTopCards();
    renderPanel(activeTab);
    const ts=new Date(statsData.updated_at);
    document.getElementById('ts').textContent=ts.toISOString().replace('T',' ').slice(0,19)+' UTC';
  }catch(e){
    document.getElementById('filter-status').textContent='統計載入失敗: '+e.message;
  }
}

function renderTopCards(){
  if(!statsData) return;
  document.getElementById('c-total').textContent=statsData.total.toLocaleString();
  const pmap={};
  statsData.purpose.forEach(p=>pmap[p.label]=p.count);
  document.getElementById('c-payload').textContent=(pmap['有效載荷']||0).toLocaleString();
  document.getElementById('c-debris') .textContent=(pmap['碎片']    ||0).toLocaleString();
  document.getElementById('c-rocket') .textContent=(pmap['火箭體']  ||0).toLocaleString();
  const emap={};
  statsData.era.forEach(e=>emap[e.label]=e.count);
  document.getElementById('c-new').textContent=(emap['< 1 年']||0).toLocaleString();
}

// ── 面板模式切換 ────────────────────────────────────────────────────────────
function setPanelMode(mode){
  panelMode=mode;
  const showTabs=(mode==='tabs');
  document.getElementById('tabs').style.display=showTabs?'flex':'none';
  document.getElementById('panel-back').style.display=showTabs?'none':'block';
  if(showTabs) renderPanel(activeTab);
}

function backToTabs(){
  // 取消 conjunction card active 狀態
  document.getElementById('conj-card').classList.remove('active-card');
  setPanelMode('tabs');
}

function switchTab(tab){
  activeTab=tab;
  if(panelMode!=='tabs') setPanelMode('tabs');
  document.querySelectorAll('.tab').forEach(t=>{
    t.classList.toggle('active',t.dataset.tab===tab);
  });
  renderPanel(tab);
}

function renderPanel(tab){
  if(!statsData||panelMode!=='tabs') return;
  const body=document.getElementById('panel-body');
  body.innerHTML='';
  const rows=statsData[tab]||[];
  const maxCount=rows.length?rows[0].count:1;
  rows.forEach(row=>{
    const isActive=(activeFtype===tab&&activeFval===row.label);
    const color=getColor(tab,row.label);
    const pct=Math.round(row.count/maxCount*100);
    const btn=document.createElement('button');
    btn.className='stat-row'+(isActive?' active':'');
    btn.innerHTML=
      '<span class="stat-dot" style="background:'+color+'"></span>'+
      '<span class="stat-label" title="'+row.label+'">'+row.label+'</span>'+
      '<span class="stat-cnt">'+row.count.toLocaleString()+'</span>';
    btn.addEventListener('click',()=>filterGlobe(tab,row.label,color));
    const barWrap=document.createElement('div');
    barWrap.style.cssText='padding:0 12px 3px 27px;width:100%';
    barWrap.innerHTML='<div class="stat-bar-wrap"><div class="stat-bar" style="width:'+pct+'%;background:'+color+'"></div></div>';
    const wrap=document.createElement('div');
    wrap.appendChild(btn); wrap.appendChild(barWrap);
    body.appendChild(wrap);
  });
}

async function filterGlobe(ftype,fval,color){
  if(activeFtype===ftype&&activeFval===fval){
    activeFtype=activeFval=null;
    satDs.entities.removeAll(); entMap.clear();
    document.getElementById('filter-status').textContent='已清除篩選';
    renderPanel(activeTab);
    return;
  }
  activeFtype=ftype; activeFval=fval;
  renderPanel(activeTab);
  document.getElementById('filter-status').textContent='載入中：'+fval+' …';
  try{
    const url='/api/positions?ftype='+encodeURIComponent(ftype)+'&fval='+encodeURIComponent(fval);
    const r=await fetch(url);
    if(!r.ok) throw new Error('HTTP '+r.status);
    const d=await r.json();
    renderEntities(d.satellites,ftype);
    const vflag=d.vectorized?'，向量化':'';
    const elapsed=d.elapsed_sec?'（'+d.elapsed_sec+' s'+vflag+'）':'';
    document.getElementById('filter-status').textContent=
      '顯示 '+d.count+' / '+d.total_matched+' 顆'+elapsed+' — '+fval;
  }catch(e){
    document.getElementById('filter-status').textContent='載入失敗: '+e.message;
  }
}

function renderEntities(sats,ftype){
  satDs.entities.removeAll(); entMap.clear();
  sats.forEach(s=>{
    const pos=Cesium.Cartesian3.fromDegrees(s.lon,s.lat,s.alt_km*1000);
    const col=Cesium.Color.fromCssColorString(s.color||'#78909C');
    const ent=satDs.entities.add({
      id:'sat_'+s.norad_id, name:s.name, position:pos,
      point:{
        pixelSize:5, color:col,
        outlineColor:Cesium.Color.WHITE.withAlpha(0.25), outlineWidth:1,
        scaleByDistance:new Cesium.NearFarScalar(5e5,1.8,1e7,0.7),
      },
      label:{
        text:s.name, font:'10px Tahoma,sans-serif',
        fillColor:Cesium.Color.WHITE, outlineColor:Cesium.Color.BLACK, outlineWidth:2,
        style:Cesium.LabelStyle.FILL_AND_OUTLINE,
        pixelOffset:new Cesium.Cartesian2(0,-14),
        distanceDisplayCondition:new Cesium.DistanceDisplayCondition(0,1.5e6),
      },
      description:new Cesium.ConstantProperty(_buildDesc(s)),
    });
    entMap.set(s.norad_id,ent);
  });
}

function _buildDesc(s){
  return '<style>body{background:#fff;color:#212121;margin:0}</style>'
    +'<div style="font-family:Tahoma,sans-serif;font-size:13px;padding:6px">'
    +'<p style="color:#1565c0;font-weight:bold;margin-bottom:6px">'+s.name+'</p>'
    +'<table style="border-collapse:collapse;color:#333;font-size:12px">'
    +'<tr><td style="padding:2px 10px 2px 0;color:#555">NORAD</td><td><b>'+s.norad_id+'</b></td></tr>'
    +'<tr><td style="padding:2px 10px 2px 0;color:#555">國家</td><td>'+(s.country||'—')+'</td></tr>'
    +'<tr><td style="padding:2px 10px 2px 0;color:#555">用途</td><td>'+(s.purpose||'—')+'</td></tr>'
    +'<tr><td style="padding:2px 10px 2px 0;color:#555">年代</td><td>'+(s.era||'—')+'</td></tr>'
    +(s.constellation&&s.constellation!=='—'?'<tr><td style="padding:2px 10px 2px 0;color:#555">星座</td><td>'+s.constellation+'</td></tr>':'')
    +'<tr><td style="padding:2px 10px 2px 0;color:#555">高度</td><td><b>'+(s.alt_km!=null?s.alt_km:'—')+' km</b></td></tr>'
    +'</table></div>';
}

// ═══════════════════════════════════════════════════════════════════════════
// 2-A  即時近距離配對
// ═══════════════════════════════════════════════════════════════════════════

async function loadConjunctions(){
  document.getElementById('c-conj').textContent='...';
  try{
    const r=await fetch('/api/conjunctions');
    if(!r.ok) throw new Error('HTTP '+r.status);
    conjData=await r.json();
    document.getElementById('c-conj').textContent=conjData.count.toLocaleString();
    if(panelMode==='conjunctions') renderConjList();
  }catch(e){
    document.getElementById('c-conj').textContent='err';
    console.warn('接近事件載入失敗',e);
  }
}

function toggleConjPanel(){
  if(panelMode==='conjunctions'){
    document.getElementById('conj-card').classList.remove('active-card');
    setPanelMode('tabs');
    return;
  }
  document.getElementById('conj-card').classList.add('active-card');
  clearSearchUI();
  setPanelMode('conjunctions');
  renderConjList();
}

function renderConjList(){
  const body=document.getElementById('panel-body');
  body.innerHTML='';
  if(!conjData){
    body.innerHTML='<div style="padding:12px;color:#484f58;font-size:11px">接近事件載入中...</div>';
    return;
  }
  if(conjData.error){
    body.innerHTML='<div style="padding:12px;color:#FF5722;font-size:11px">'+conjData.error+'</div>';
    return;
  }
  if(!conjData.pairs||!conjData.pairs.length){
    body.innerHTML='<div style="padding:12px;color:#484f58;font-size:11px">'
      +'目前 '+conjData.threshold_km+' km 閾值內無配對（掃描 '+(conjData.total_scanned||0)+' 顆）</div>';
    return;
  }
  // 注釋
  const note=document.createElement('div');
  note.className='conj-note';
  note.textContent='閾值 '+conjData.threshold_km+' km | 掃描 '+conjData.total_scanned
    +' 顆 | 耗時 '+conjData.elapsed_sec+' s'+(conjData.vectorized?' | ⚡向量化':'');
  body.appendChild(note);

  conjData.pairs.forEach(p=>{
    const row=document.createElement('div');
    row.className='conj-row';
    const purposeColor={有效載荷:'#4CAF50',碎片:'#FF9800',火箭體:'#9E9E9E'}[p.primary_purpose]||'#78909C';
    const sPurposeColor={有效載荷:'#4CAF50',碎片:'#FF9800',火箭體:'#9E9E9E'}[p.secondary_purpose]||'#78909C';
    row.innerHTML=
      '<div style="display:flex;align-items:center;gap:6px">'
      +'<span class="conj-km">'+p.miss_km.toFixed(2)+' km</span>'
      +'<span class="conj-names" title="'+p.primary_name+' + '+p.secondary_name+'">'
      +'<span style="color:'+purposeColor+'">●</span> '+p.primary_name
      +' <span style="color:#484f58">+</span>'
      +' <span style="color:'+sPurposeColor+'">●</span> '+p.secondary_name
      +'</span></div>'
      +'<div style="font-size:9px;color:#484f58;margin-top:2px">'
      +p.primary_alt_km+' km / '+p.secondary_alt_km+' km</div>';
    row.addEventListener('click',()=>flyToPair(p));
    body.appendChild(row);
  });

  const footer=document.createElement('div');
  footer.className='conj-note';
  footer.innerHTML='<span style="color:#484f58">詳細接近分析（TCA/Pc）→ conjunction_pipeline.py</span>';
  body.appendChild(footer);
}

function flyToPair(pair){
  satDs.entities.removeAll(); entMap.clear();
  const posA=Cesium.Cartesian3.fromDegrees(pair.primary_lon,pair.primary_lat,pair.primary_alt_km*1000);
  const posB=Cesium.Cartesian3.fromDegrees(pair.secondary_lon,pair.secondary_lat,pair.secondary_alt_km*1000);
  const colA=Cesium.Color.fromCssColorString('#FF6B6B');
  const colB=Cesium.Color.fromCssColorString('#FFB86C');
  satDs.entities.add({id:'ca_a',name:pair.primary_name,position:posA,
    point:{pixelSize:10,color:colA,outlineColor:Cesium.Color.WHITE,outlineWidth:2},
    description:new Cesium.ConstantProperty(
      '<div style="font-family:Tahoma;padding:6px;font-size:13px;color:#212121">'
      +'<b>'+pair.primary_name+'</b><br>NORAD '+pair.primary_norad
      +'<br>高度 '+pair.primary_alt_km+' km</div>')
  });
  satDs.entities.add({id:'ca_b',name:pair.secondary_name,position:posB,
    point:{pixelSize:10,color:colB,outlineColor:Cesium.Color.WHITE,outlineWidth:2},
    description:new Cesium.ConstantProperty(
      '<div style="font-family:Tahoma;padding:6px;font-size:13px;color:#212121">'
      +'<b>'+pair.secondary_name+'</b><br>NORAD '+pair.secondary_norad
      +'<br>高度 '+pair.secondary_alt_km+' km</div>')
  });
  satDs.entities.add({polyline:{
    positions:[posA,posB], width:2,
    material:new Cesium.ColorMaterialProperty(Cesium.Color.fromCssColorString('#FF6B6B').withAlpha(0.7)),
    arcType:Cesium.ArcType.NONE,
  }});
  const midLat=(pair.primary_lat+pair.secondary_lat)/2;
  const midLon=(pair.primary_lon+pair.secondary_lon)/2;
  const midAlt=(pair.primary_alt_km+pair.secondary_alt_km)/2*1000;
  viewer.camera.flyTo({
    destination:Cesium.Cartesian3.fromDegrees(midLon,midLat,midAlt+800000),
    duration:2,
  });
  document.getElementById('filter-status').textContent=
    pair.primary_name+' + '+pair.secondary_name+' — '+pair.miss_km.toFixed(2)+' km';
}

// ═══════════════════════════════════════════════════════════════════════════
// 4-A  物件搜尋
// ═══════════════════════════════════════════════════════════════════════════

document.getElementById('search-input').addEventListener('input', function(){
  clearTimeout(_searchTimer);
  const q=this.value.trim();
  document.getElementById('search-clear').style.display=q?'':'none';
  if(q.length<2){ clearSearchUI(); return; }
  _searchTimer=setTimeout(()=>doSearch(q),400);
});

async function doSearch(q){
  try{
    const r=await fetch('/api/search?q='+encodeURIComponent(q));
    if(!r.ok) return;
    const d=await r.json();
    showSearchResults(d.results);
  }catch(e){ console.warn('搜尋失敗',e); }
}

function showSearchResults(results){
  document.getElementById('conj-card').classList.remove('active-card');
  setPanelMode('search');
  const body=document.getElementById('panel-body');
  body.innerHTML='';
  if(!results||!results.length){
    body.innerHTML='<div style="padding:12px;color:#484f58;font-size:11px">未找到符合結果</div>';
    return;
  }
  results.forEach(s=>{
    const color=getColor('country',s.country);
    const row=document.createElement('button');
    row.className='stat-row';
    row.innerHTML=
      '<span class="stat-dot" style="background:'+color+'"></span>'
      +'<span class="stat-label" title="'+s.name+'">'+s.name+'</span>'
      +'<span style="font-size:10px;color:#484f58;min-width:54px;text-align:right">'+s.norad_id+'</span>';
    row.addEventListener('click',()=>flyToSat(s));
    body.appendChild(row);
  });
}

function clearSearch(){
  document.getElementById('search-input').value='';
  clearSearchUI();
}

function clearSearchUI(){
  document.getElementById('search-clear').style.display='none';
  clearTimeout(_searchTimer);
  if(panelMode==='search') setPanelMode('tabs');
}

function flyToSat(s){
  if(s.lat==null) return;
  satDs.entities.removeAll(); entMap.clear();
  const pos=Cesium.Cartesian3.fromDegrees(s.lon,s.lat,s.alt_km*1000);
  const col=Cesium.Color.fromCssColorString(getColor('country',s.country));
  const ent=satDs.entities.add({
    id:'sat_'+s.norad_id, name:s.name, position:pos,
    point:{pixelSize:10,color:col,outlineColor:Cesium.Color.WHITE,outlineWidth:2},
    description:new Cesium.ConstantProperty(_buildDesc(s)),
  });
  viewer.flyTo(ent,{duration:2});
  viewer.selectedEntity=ent;
  document.getElementById('filter-status').textContent=s.name+' (#'+s.norad_id+')';
  activeFtype=activeFval=null;
}

// ── 向量圖層管理 ──────────────────────────────────────────────────────────
async function toggleLayer(type,cb){
  if(type==='borders'){
    if(borderDs){ borderDs.show=cb.checked; return; }
    if(cb.checked) await loadBordersLayer();
    return;
  }
  const ref={ssn:{get:()=>ssnDs,set:v=>{ssnDs=v;},url:'/api/layers/ssn_stations'}}[type];
  if(!ref) return;
  if(cb.checked){
    try{
      const ds=await Cesium.GeoJsonDataSource.load(ref.url);
      [...ds.entities.values].forEach(ent=>{
        const props=ent.properties;
        const stationType=props.type&&props.type.getValue()||'雷達';
        const status=props.status&&props.status.getValue()||'active';
        const hexCol=SSN_TYPE_COLORS[stationType]||'#00E5FF';
        const col=Cesium.Color.fromCssColorString(hexCol);
        const isRetired=(status==='decommissioned');
        ent.billboard=undefined;
        ent.point=new Cesium.PointGraphics({
          pixelSize:isRetired?5:9,
          color:isRetired?col.withAlpha(0.45):col,
          outlineColor:Cesium.Color.BLACK.withAlpha(0.7),outlineWidth:1,
        });
        const nameVal=props.name&&props.name.getValue()||'';
        ent.label=new Cesium.LabelGraphics({
          text:nameVal,font:'10px Tahoma,sans-serif',
          fillColor:isRetired?Cesium.Color.fromCssColorString('#9E9E9E'):Cesium.Color.WHITE,
          outlineColor:Cesium.Color.BLACK,outlineWidth:2,
          style:Cesium.LabelStyle.FILL_AND_OUTLINE,
          pixelOffset:new Cesium.Cartesian2(0,-13),
          distanceDisplayCondition:new Cesium.DistanceDisplayCondition(0,isRetired?2e6:7e6),
          show:true,
        });
        const locVal=props.location&&props.location.getValue()||'';
        const notesVal=props.notes&&props.notes.getValue()||'';
        ent.description=new Cesium.ConstantProperty(
          '<div style="font-family:Tahoma,sans-serif;font-size:13px;padding:8px;background:#1e1e1e;color:#e8e8e8">'
          +'<p style="color:#8ab4f8;font-weight:bold;font-size:14px;margin:0 0 8px">'+nameVal+'</p>'
          +'<table style="border-collapse:collapse;width:100%;font-size:12px">'
          +'<tr><td style="color:#aaa;padding:2px 10px 2px 0">類型</td><td style="color:#e8e8e8"><b>'+stationType+'</b></td></tr>'
          +'<tr><td style="color:#aaa;padding:2px 10px 2px 0">位置</td><td style="color:#e8e8e8">'+locVal+'</td></tr>'
          +'<tr><td style="color:#aaa;padding:2px 10px 2px 0">狀態</td><td>'
          +(isRetired?'<span style="color:#f44336">已除役</span>':'<span style="color:#4caf50">運作中</span>')
          +'</td></tr>'+(notesVal?'<tr><td style="color:#aaa;padding:2px 10px 2px 0;vertical-align:top">備註</td><td style="color:#ccc">'+notesVal+'</td></tr>':'')
          +'</table></div>'
        );
        const posCart=ent.position.getValue(Cesium.JulianDate.now());
        if(posCart){
          const carto=Cesium.Cartographic.fromCartesian(posCart);
          const lon=Cesium.Math.toDegrees(carto.longitude);
          const lat=Cesium.Math.toDegrees(carto.latitude);
          const h=20000;
          ent.position=new Cesium.ConstantPositionProperty(Cesium.Cartesian3.fromDegrees(lon,lat,h));
          ds.entities.add({polyline:{
            positions:Cesium.Cartesian3.fromDegreesArrayHeights([lon,lat,0,lon,lat,h]),
            width:1,
            material:new Cesium.ColorMaterialProperty(
              (isRetired?Cesium.Color.fromCssColorString('#9E9E9E'):col).withAlpha(isRetired?0.30:0.55)),
            clampToGround:false,arcType:Cesium.ArcType.NONE,
            distanceDisplayCondition:new Cesium.DistanceDisplayCondition(0,7e6),
          }});
        }
      });
      await viewer.dataSources.add(ds);
      ref.set(ds);
    }catch(e){
      cb.checked=false;
      console.warn('向量圖層載入失敗: '+type,e);
    }
  }else{
    const ds=ref.get();
    if(ds) viewer.dataSources.remove(ds,true);
    ref.set(null);
  }
}

// ── 啟動 ──────────────────────────────────────────────────────────────────
async function init(){
  await initCesium();
  await loadBordersLayer();
  await loadStats();
  loadConjunctions();   // 背景載入，不 await
}

init().catch(e=>{
  document.getElementById('loading').textContent='初始化失敗: '+e.message;
  console.error(e);
});

window.switchTab=switchTab;
window.toggleLayer=toggleLayer;
window.toggleConjPanel=toggleConjPanel;
window.clearSearch=clearSearch;
window.backToTabs=backToTabs;

} // end startApp
</script>
</body>
</html>"""


# ─────────────────────────────────────────────────────────────────────────────
# Flask 應用程式
# ─────────────────────────────────────────────────────────────────────────────
import json as _json


def create_app() -> Flask:
    app = Flask(__name__)
    CORS(app, resources={r"/api/*": {"origins": "*"}})

    @app.get("/")
    def index():
        return _HTML, 200, {"Content-Type": "text/html; charset=utf-8"}

    # ── 本機 Cesium ──────────────────────────────────────────────────────────

    @app.get("/cesium/<path:filename>")
    def cesium_static(filename: str):
        safe = (_CESIUM_LOCAL_DIR / filename).resolve()
        if not str(safe).startswith(str(_CESIUM_LOCAL_DIR.resolve())):
            return make_response("Forbidden", 403)
        if not safe.is_file():
            return make_response(f"Cesium asset not found: {filename}", 404)
        return send_from_directory(str(_CESIUM_LOCAL_DIR), filename)

    # ── Globe 貼圖 ───────────────────────────────────────────────────────────

    @app.get("/api/globe_texture")
    def api_globe_texture():
        data = _get_globe_texture()
        if data is None:
            return make_response("Globe 貼圖無法取得", 503)
        ct = "image/png" if data[:4] == b"\x89PNG" else "image/jpeg"
        resp = make_response(data)
        resp.headers["Content-Type"]  = ct
        resp.headers["Cache-Control"] = "public, max-age=604800"
        return resp

    # ── 統計 ─────────────────────────────────────────────────────────────────

    @app.get("/api/stats")
    def api_stats():
        return jsonify(get_stats())

    # ── 1-A  向量化 SGP4 全量傳播 ────────────────────────────────────────────

    @app.get("/api/positions")
    def api_positions():
        ftype = request.args.get("ftype", "").strip()
        fval  = request.args.get("fval",  "").strip()
        VALID = {"country", "purpose", "era", "constellation"}
        if ftype not in VALID or not fval:
            return jsonify({"error": "ftype 必須為 country/purpose/era/constellation，且 fval 不可空白"}), 400

        idx = get_sat_index()
        if ftype == "constellation":
            matched = [n for n, i in idx.items() if i.get("constellation") == fval]
        else:
            matched = [n for n, i in idx.items() if i.get(ftype) == fval]

        total = len(matched)
        t0    = time.monotonic()
        logger.info("向量化傳播 %d 顆（%s=%s）", total, ftype, fval)

        positions = propagate_batch(matched, idx)
        elapsed   = time.monotonic() - t0
        logger.info("傳播完成 %d 顆，耗時 %.2f s（%s）",
                    total, elapsed, "SatrecArray" if _HAS_SATREC_ARRAY else "sequential")

        color   = get_color(ftype, fval)
        results = []
        for nid, pos in zip(matched, positions):
            if pos is None:
                continue
            lat, lon, alt = pos
            info = idx[nid]
            results.append({
                "norad_id":      nid,
                "name":          info["name"],
                "country":       info["country"],
                "purpose":       info["purpose"],
                "era":           info["era"],
                "constellation": info["constellation"] or "—",
                "color":         color,
                "lat":           round(lat, 4),
                "lon":           round(lon, 4),
                "alt_km":        round(alt, 1),
            })

        return jsonify({
            "ftype":         ftype,
            "fval":          fval,
            "count":         len(results),
            "total_matched": total,
            "sampled":       False,
            "elapsed_sec":   round(elapsed, 3),
            "vectorized":    _HAS_SATREC_ARRAY,
            "satellites":    results,
            "timestamp":     datetime.now(timezone.utc).isoformat(),
        })

    # ── 單顆衛星查詢（供搜尋結果使用）──────────────────────────────────────

    @app.get("/api/position/<int:norad_id>")
    def api_position_single(norad_id: int):
        idx = get_sat_index()
        info = idx.get(norad_id)
        if info is None:
            return jsonify({"error": f"NORAD {norad_id} 不在索引中"}), 404
        pos = propagate_now(info["line1"], info["line2"])
        result: dict[str, Any] = {
            "norad_id":      norad_id,
            "name":          info["name"],
            "country":       info["country"],
            "purpose":       info["purpose"],
            "era":           info["era"],
            "constellation": info["constellation"] or "—",
        }
        if pos:
            result["lat"]    = round(pos[0], 4)
            result["lon"]    = round(pos[1], 4)
            result["alt_km"] = round(pos[2], 1)
        return jsonify(result)

    # ── 2-A  近距離配對掃描 ──────────────────────────────────────────────────

    @app.get("/api/conjunctions")
    def api_conjunctions():
        try:
            threshold = float(request.args.get("threshold_km", _CONJ_THRESHOLD_KM))
            threshold = max(1.0, min(threshold, 500.0))
        except ValueError:
            threshold = _CONJ_THRESHOLD_KM
        try:
            max_pairs = int(request.args.get("max_pairs", 200))
            max_pairs = max(10, min(max_pairs, 2000))
        except ValueError:
            max_pairs = 200

        global _conj_cache, _conj_loaded_at
        if (_conj_cache is None
                or (time.monotonic() - _conj_loaded_at) > _CONJ_TTL
                or _conj_cache.get("threshold_km") != threshold):
            _conj_cache     = build_conjunction_summary(threshold, max_pairs)
            _conj_loaded_at = time.monotonic()

        resp = make_response(_json.dumps(_conj_cache, ensure_ascii=False).encode("utf-8"))
        resp.headers["Content-Type"]  = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = f"public, max-age={_CONJ_TTL}"
        return resp

    # ── 4-A  物件搜尋 ────────────────────────────────────────────────────────

    @app.get("/api/search")
    def api_search():
        q = request.args.get("q", "").strip()
        if len(q) < 2:
            return jsonify({"results": [], "count": 0, "query": q})

        idx  = get_sat_index()
        q_up = q.upper()
        matches: list[dict[str, Any]] = []

        # 精確 NORAD ID 比對（優先）
        if q.isdigit():
            nid = int(q)
            if nid in idx:
                matches.append({"norad_id": nid, **idx[nid], "score": 0})

        # 名稱子串比對
        for nid, info in idx.items():
            if q_up in info["name"].upper():
                if not any(m["norad_id"] == nid for m in matches):
                    matches.append({"norad_id": nid, **info, "score": 1})
            if len(matches) >= 60:
                break

        matches.sort(key=lambda x: (x["score"], x["name"]))
        top = matches[:20]

        # 批次傳播 top-20 位置
        nids      = [m["norad_id"] for m in top]
        positions = propagate_batch(nids, idx)

        results = []
        for m, pos in zip(top, positions):
            r: dict[str, Any] = {
                "norad_id":      m["norad_id"],
                "name":          m["name"],
                "country":       m["country"],
                "purpose":       m["purpose"],
                "era":           m["era"],
                "constellation": m["constellation"] or "—",
            }
            if pos:
                r["lat"]    = round(pos[0], 4)
                r["lon"]    = round(pos[1], 4)
                r["alt_km"] = round(pos[2], 1)
            results.append(r)

        return jsonify({"results": results, "count": len(results), "query": q})

    # ── 向量圖層 ─────────────────────────────────────────────────────────────

    @app.get("/api/layers/borders")
    def api_layer_borders():
        global _borders_cache
        if _borders_cache is None:
            if _BORDERS_LOCAL.exists():
                try:
                    _borders_cache = _BORDERS_LOCAL.read_bytes()
                    logger.info("國界 GeoJSON 本地: %d bytes", len(_borders_cache))
                except Exception as exc:
                    logger.warning("國界讀取失敗: %s", exc)
            if _borders_cache is None:
                try:
                    logger.info("下載 Natural Earth 國界: %s", _NE_BORDERS_URL)
                    r = requests.get(_NE_BORDERS_URL, timeout=30,
                                     headers={"User-Agent": "ATRDC-TLE-Tracker/1.0"})
                    r.raise_for_status()
                    raw = _json.loads(r.content)
                    for feat in raw.get("features", []):
                        feat["properties"] = {}
                    _borders_cache = _json.dumps(raw, separators=(",", ":")).encode("utf-8")
                    try:
                        _BORDERS_LOCAL.parent.mkdir(parents=True, exist_ok=True)
                        _BORDERS_LOCAL.write_bytes(_borders_cache)
                    except Exception:
                        pass
                except Exception as exc:
                    logger.warning("國界下載失敗: %s", exc)
                    _borders_cache = b'{"type":"FeatureCollection","features":[]}'
        resp = make_response(_borders_cache)
        resp.headers["Content-Type"]  = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = "public, max-age=86400"
        return resp

    @app.get("/api/layers/ssn_stations")
    def api_layer_ssn_stations():
        resp = make_response(
            _json.dumps(_SSN_STATIONS_GEOJSON, ensure_ascii=False).encode("utf-8"))
        resp.headers["Content-Type"]  = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = "public, max-age=3600"
        return resp

    @app.errorhandler(Exception)
    def handle_error(err):
        from werkzeug.exceptions import HTTPException
        if isinstance(err, HTTPException):
            return err
        logger.exception("未預期例外")
        return jsonify({"error": str(err)}), 500

    return app


app = create_app()

if __name__ == "__main__":
    logger.info(
        "Scenario 04-advanced01 (向量化 SGP4 + 近距離掃描 + 搜尋) 啟動 — http://%s:%d",
        HOST, PORT,
    )
    logger.info(
        "SatrecArray=%s  KD-tree=%s  接近閾值=%.0f km  快取 TTL=%d s",
        _HAS_SATREC_ARRAY, _HAS_KDTREE, _CONJ_THRESHOLD_KM, _CONJ_TTL,
    )
    logger.info("預熱衛星索引…")
    get_sat_index()
    get_stats()
    app.run(host=HOST, port=PORT, debug=True)
