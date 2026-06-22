#!/usr/bin/env python3
"""
scenario04-Cesium-advanced02.py — 太空態勢儀表板（進階版）+ 台北覆蓋分析 + 時間軸
====================================================================================
Port 5013

以 scenario04-Cesium-advanced01-2.py 為基礎，台北覆蓋頁面新增時間軸功能：
  - 可回溯過去 30 天（使用當時最新 TLE 從 DB 查詢）
  - 可預測未來 30 天（使用現行 TLE + SGP4 外推）
  - 時間軸滑桿（步長 15 分鐘）+ 「現在」按鈕
  - 歷史 / 現在 / 預測 模式指示
  - 處於「現在」時每 60 秒自動更新

新增路由：
  GET /taipei                  → Cesium ION 2D 台北覆蓋頁面
  GET /api/taipei_coverage     → 當前各類衛星位置 + 仰角方位
  GET /api/taipei_passes       → 24 小時過頂預報（四類）

四類衛星：
  1. US_EO   — 美國商用光學（Vantor/Maxar WorldView、Planet SkySat/Pelican）
  2. CN_COMM — 中國商用光學（SuperView、高分、吉林）
  3. CN_MIL  — 中國軍用偵察（遙感 Yaogan）
  4. TW_TASA — 台灣 TASA（Formosat-5/-7/-8、COSMIC-2）

Cesium 2D 穩定關鍵（避免畫面向下飄移）：
  - 直接以 SCENE2D 初始化，不走 morph
  - camera.setView()（同步），不用 flyTo()（非同步動畫）
  - viewer.clock.shouldAnimate = false
  - 停用 screenSpaceCameraController 的 enableRotate / enableTilt
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
DEFAULT_DB  = str(_SCRIPT_DIR / "space_db.duckdb")
RAW_TABLE   = "raw_tle_archive"
META_TABLE  = "sat_n2yo_metadata"

load_dotenv(_SCRIPT_DIR / ".env")

logging.basicConfig(
    level=os.getenv("LOG_LEVEL", "INFO").upper(),
    format="%(asctime)s %(levelname)s %(name)s %(message)s",
)
logger = logging.getLogger("scenario04_cesium_adv01_2")

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
PORT                  = int(os.getenv("PORT", "5013"))
STATS_TTL             = int(os.getenv("STATS_TTL", "600"))
_CONJ_THRESHOLD_KM    = float(os.getenv("CONJ_THRESHOLD_KM", "10.0"))
_CONJ_TTL             = int(os.getenv("CONJ_TTL", "120"))

# ── 台北覆蓋分析常數 ──────────────────────────────────────────────────────────
_TAIPEI_LAT      = 25.0330
_TAIPEI_LON      = 121.5654
_TAIPEI_H_KM     = 0.01
_COVER_KM        = 2000.0
_MASK_DEG        = 5.0
CESIUM_ION_TOKEN = os.getenv("CESIUM_ION_TOKEN", "")

_OVERPASS_CATS: dict[str, dict] = {
    "US_EO": {
        "label":    "美國商用光學衛星",
        "sublabel": "Vantor/Maxar · Planet SkySat/Pelican",
        "color":    "#4488FF",
        "kw":       ["WORLDVIEW", "GEOEYE", "LEGION", "SKYSAT", "PELICAN"],
    },
    "CN_COMM": {
        "label":    "中國商用光學衛星",
        "sublabel": "SuperView · 高分 · 吉林",
        "color":    "#FF9800",
        "kw":       ["SUPERVIEW", "JILIN", "ZHUHAI", "GAOFEN"],
    },
    "CN_MIL": {
        "label":    "中國軍用偵察衛星",
        "sublabel": "遙感 Yaogan",
        "color":    "#F44336",
        "kw":       ["YAOGAN", "JIANBING"],
    },
    "TW_TASA": {
        "label":    "台灣 TASA 衛星",
        "sublabel": "Formosat-5 / -7 / -8",
        "color":    "#00E5FF",
        "kw":       [
            "FORMOSAT-5", "FORMOSAT 5", "FORMOSAT5",
            "FORMOSAT-7", "FORMOSAT 7", "FORMOSAT7",
            "FORMOSAT-8", "FORMOSAT 8", "FORMOSAT8",
            "COSMIC-2", "COSMIC2",
        ],
    },
}

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

_db_info_cache:     dict[str, Any] = {}
_db_info_loaded_at: float = 0.0
_DB_INFO_TTL = 300


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
    alt_name = "space_db_slim.duckdb" if DB_PATH.name == "space_db.duckdb" else "space_db.duckdb"
    alt = DB_PATH.parent / alt_name
    if alt.exists():
        logger.warning("DB %s 不存在，改用 %s", DB_PATH.name, alt.name)
        return alt
    logger.error("找不到資料庫: %s", DB_PATH)
    return None


def _tle_select_sql(con: duckdb.DuckDBPyConnection, extra_where: str = "") -> str:
    """動態偵測 raw_tle_archive 欄位，生成相容 slim/full DB 的查詢 SQL。"""
    actual = {r[0] for r in con.execute(f"DESCRIBE {RAW_TABLE}").fetchall()}
    obj_expr = (
        "COALESCE(r.object_name, 'NORAD-' || CAST(r.norad_id AS VARCHAR))"
        if "object_name" in actual
        else "'NORAD-' || CAST(r.norad_id AS VARCHAR)"
    )
    has_lines = "line1" in actual and "line2" in actual
    line_sel   = "r.line1, r.line2," if has_lines else "NULL AS line1, NULL AS line2,"
    line_where = "r.line1 IS NOT NULL AND r.line2 IS NOT NULL" if has_lines else "1=1"

    parts = [line_where]
    if extra_where:
        parts.append(extra_where)

    return f"""
        SELECT
            r.norad_id,
            {obj_expr} AS object_name,
            {line_sel}
            m.source_code, m.launch_date, m.intl_code
        FROM {RAW_TABLE} r
        LEFT JOIN {META_TABLE} m ON r.norad_id = m.norad_id
        WHERE {" AND ".join(parts)}
        QUALIFY ROW_NUMBER() OVER (
            PARTITION BY r.norad_id ORDER BY r.epoch_utc DESC
        ) = 1
    """


def build_sat_index() -> dict[int, dict[str, Any]]:
    db = _resolve_db()
    if db is None:
        return {}
    logger.info("建立衛星索引中…")
    t0 = time.monotonic()
    try:
        with duckdb.connect(str(db), read_only=True) as con:
            rows = con.execute(_tle_select_sql(con)).fetchall()
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


def get_db_info() -> dict[str, Any]:
    global _db_info_cache, _db_info_loaded_at
    if _db_info_cache and (time.monotonic() - _db_info_loaded_at) < _DB_INFO_TTL:
        return _db_info_cache

    db = _resolve_db()
    if db is None:
        return {"error": "資料庫不存在"}

    try:
        mtime_ts = db.stat().st_mtime
        db_updated_at = datetime.fromtimestamp(mtime_ts, tz=timezone.utc).isoformat()
        db_size_mb = round(db.stat().st_size / 1024**2, 1)

        with duckdb.connect(str(db), read_only=True) as con:
            actual_cols = {r[0] for r in con.execute("DESCRIBE raw_tle_archive").fetchall()}
            has_lines = "line1" in actual_cols and "line2" in actual_cols

            where = "WHERE line1 IS NOT NULL AND line2 IS NOT NULL" if has_lines else ""
            row = con.execute(f"""
                SELECT
                    COUNT(*)                 AS total_records,
                    COUNT(DISTINCT norad_id) AS valid_sat_count,
                    MIN(epoch_utc)           AS epoch_min,
                    MAX(epoch_utc)           AS epoch_max
                FROM raw_tle_archive {where}
            """).fetchone()

        def _iso(v: Any) -> str | None:
            if v is None:
                return None
            if hasattr(v, "isoformat"):
                return v.isoformat()
            return str(v)

        result: dict[str, Any] = {
            "db_name":         db.name,
            "db_size_mb":      db_size_mb,
            "db_updated_at":   db_updated_at,
            "has_tle_lines":   has_lines,
            "total_records":   int(row[0]) if row and row[0] else 0,
            "valid_sat_count": int(row[1]) if row and row[1] else 0,
            "epoch_min":       _iso(row[2]) if row else None,
            "epoch_max":       _iso(row[3]) if row else None,
        }
    except Exception as exc:
        result = {"error": str(exc), "db_name": db.name}

    _db_info_cache = result
    _db_info_loaded_at = time.monotonic()
    return result


# ─────────────────────────────────────────────────────────────────────────────
# 1-A  向量化 SGP4
# ─────────────────────────────────────────────────────────────────────────────

def eci_to_llh_batch(r_arr: np.ndarray, t: datetime) -> np.ndarray:
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
    n   = len(nids)
    jd0, fr0 = jday(t.year, t.month, t.day, t.hour, t.minute,
                    t.second + t.microsecond * 1e-6)

    if _HAS_SATREC_ARRAY and n > 1:
        try:
            sats = _SatrecArray([Satrec.twoline2rv(l1, l2) for l1, l2 in zip(line1s, line2s)])
            e_raw, r_raw, _ = sats.sgp4(np.array([jd0]), np.array([fr0]))
            return e_raw[:, 0].astype(int), r_raw[:, 0, :]
        except Exception as exc:
            logger.debug("SatrecArray 傳播失敗，退回逐顆: %s", exc)

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
# 台北覆蓋分析：仰角計算函式
# ─────────────────────────────────────────────────────────────────────────────

def _gmst_rad(jd: float, fr: float) -> float:
    T = ((jd - 2451545.0) + fr) / 36525.0
    return np.deg2rad(
        (280.46061837 + 360.98564736629 * (jd - 2451545.0 + fr)
         + 0.000387933 * T ** 2) % 360.0)


def _observer_ecef(lat_deg: float, lon_deg: float, h_km: float) -> dict:
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    sl = float(np.sin(lat)); cl = float(np.cos(lat))
    so = float(np.sin(lon)); co = float(np.cos(lon))
    N  = R_EARTH_KM / np.sqrt(1.0 - E2 * sl ** 2)
    x0 = (N + h_km) * cl * co
    y0 = (N + h_km) * cl * so
    z0 = (N * (1.0 - E2) + h_km) * sl
    return {"x0": float(x0), "y0": float(y0), "z0": float(z0),
            "sl": sl, "cl": cl, "so": so, "co": co}


def _eci_to_elaz(
    r_eci: np.ndarray,   # (N, 3) ECI km
    jd: float, fr: float,
    obs: dict,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns (el_deg, az_deg, range_km) arrays of shape (N,)."""
    gmst = _gmst_rad(jd, fr)
    cg = np.cos(gmst); sg = np.sin(gmst)
    xe = cg * r_eci[:, 0] + sg * r_eci[:, 1]
    ye = -sg * r_eci[:, 0] + cg * r_eci[:, 1]
    ze = r_eci[:, 2]
    dx = xe - obs["x0"]; dy = ye - obs["y0"]; dz = ze - obs["z0"]
    sl = obs["sl"]; cl = obs["cl"]; so = obs["so"]; co = obs["co"]
    E_enu = -so * dx + co * dy
    N_enu = -sl * co * dx - sl * so * dy + cl * dz
    U_enu =  cl * co * dx + cl * so * dy + sl * dz
    rng   = np.sqrt(E_enu ** 2 + N_enu ** 2 + U_enu ** 2)
    safe  = np.where(rng > 0.001, rng, 0.001)
    el    = np.rad2deg(np.arcsin(np.clip(U_enu / safe, -1.0, 1.0)))
    az    = np.rad2deg(np.arctan2(E_enu, N_enu)) % 360.0
    return el, az, rng


# Pre-compute Taipei observer ECEF (used in route handlers)
_TAIPEI_OBS = _observer_ecef(_TAIPEI_LAT, _TAIPEI_LON, _TAIPEI_H_KM)


# ─────────────────────────────────────────────────────────────────────────────
# 台北覆蓋分析：衛星篩選與計算
# ─────────────────────────────────────────────────────────────────────────────

def _get_overpass_candidates(idx: dict) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {cat: [] for cat in _OVERPASS_CATS}
    for nid, info in idx.items():
        name_up = info["name"].upper()
        for cat, cfg in _OVERPASS_CATS.items():
            if any(kw in name_up for kw in cfg["kw"]):
                result[cat].append(nid)
                break
    return result


def compute_taipei_coverage(mask_deg: float = _MASK_DEG) -> dict[str, Any]:
    idx = get_sat_index()
    candidates = _get_overpass_candidates(idx)
    t = datetime.now(timezone.utc)
    jd0, fr0 = jday(t.year, t.month, t.day, t.hour, t.minute,
                    t.second + t.microsecond * 1e-6)
    obs = _TAIPEI_OBS

    all_nids = [nid for nids in candidates.values() for nid in nids]
    if not all_nids:
        return {
            "categories": {
                cat: {**cfg, "count": 0, "visible_count": 0, "satellites": []}
                for cat, cfg in _OVERPASS_CATS.items()
            },
            "timestamp": t.isoformat(),
            "mask_deg":  mask_deg,
        }

    line1s = [idx[n]["line1"] for n in all_nids]
    line2s = [idx[n]["line2"] for n in all_nids]
    err_arr, r_arr = _sgp4_propagate_raw(all_nids, line1s, line2s, t)

    valid = err_arr == 0
    el_all  = np.full(len(all_nids), -90.0)
    az_all  = np.zeros(len(all_nids))
    rng_all = np.zeros(len(all_nids))
    llh_all = np.full((len(all_nids), 3), np.nan)

    if valid.any():
        el_v, az_v, rng_v = _eci_to_elaz(r_arr[valid], jd0, fr0, obs)
        el_all[valid]  = el_v
        az_all[valid]  = az_v
        rng_all[valid] = rng_v
        llh_all[valid] = eci_to_llh_batch(r_arr[valid], t)

    nid_to_i = {nid: i for i, nid in enumerate(all_nids)}
    result_cats: dict[str, Any] = {}

    for cat, cfg in _OVERPASS_CATS.items():
        sats = []
        for nid in candidates[cat]:
            i = nid_to_i[nid]
            if err_arr[i] != 0 or np.isnan(llh_all[i, 0]):
                continue
            lat, lon, alt = llh_all[i, 0], llh_all[i, 1], llh_all[i, 2]
            if not (-500.0 < float(alt) < 80_000.0):
                continue
            el  = float(el_all[i])
            az  = float(az_all[i])
            rng = float(rng_all[i])
            sats.append({
                "norad_id":  nid,
                "name":      idx[nid]["name"],
                "lat":       round(float(lat), 4),
                "lon":       round(float(lon), 4),
                "alt_km":    round(float(alt), 1),
                "el_deg":    round(el, 2),
                "az_deg":    round(az, 2),
                "range_km":  round(rng, 1),
                "visible":   el >= mask_deg,
                "color":     cfg["color"],
            })
        result_cats[cat] = {
            "label":         cfg["label"],
            "sublabel":      cfg["sublabel"],
            "color":         cfg["color"],
            "count":         len(sats),
            "visible_count": sum(1 for s in sats if s["visible"]),
            "satellites":    sats,
        }

    return {"categories": result_cats, "timestamp": t.isoformat(), "mask_deg": mask_deg}


def predict_taipei_passes(
    hours:       float = 24.0,
    step_sec:    float = 60.0,
    mask_deg:    float = _MASK_DEG,
    max_per_cat: int   = 20,
) -> dict[str, Any]:
    idx = get_sat_index()
    candidates = _get_overpass_candidates(idx)
    t0 = datetime.now(timezone.utc)
    obs = _TAIPEI_OBS

    n_steps = int(hours * 3600 / step_sec)
    times   = [t0 + timedelta(seconds=i * step_sec) for i in range(n_steps)]

    # Build JD/FR arrays once
    jd_fr = np.array([
        jday(tt.year, tt.month, tt.day, tt.hour, tt.minute,
             tt.second + tt.microsecond * 1e-6)
        for tt in times
    ])
    jds = np.ascontiguousarray(jd_fr[:, 0])
    frs = np.ascontiguousarray(jd_fr[:, 1])

    # Pre-compute GMST for all time steps
    T_cent   = ((jds - 2451545.0) + frs) / 36525.0
    gmst_all = np.deg2rad(
        (280.46061837 + 360.98564736629 * (jds - 2451545.0 + frs)
         + 0.000387933 * T_cent ** 2) % 360.0)
    cg_all = np.cos(gmst_all)
    sg_all = np.sin(gmst_all)

    all_passes: dict[str, list] = {cat: [] for cat in _OVERPASS_CATS}

    for cat, cfg in _OVERPASS_CATS.items():
        nids = candidates[cat]
        if not nids:
            continue

        line1s = [idx[n]["line1"] for n in nids]
        line2s = [idx[n]["line2"] for n in nids]

        if not _HAS_SATREC_ARRAY:
            logger.warning("SatrecArray 不可用，跳過過頂預報 cat=%s", cat)
            continue

        try:
            sa = _SatrecArray([Satrec.twoline2rv(l1, l2) for l1, l2 in zip(line1s, line2s)])
            e_raw, r_raw, _ = sa.sgp4(jds, frs)
            # e_raw: (N_sats, N_times)  r_raw: (N_sats, N_times, 3)
        except Exception as exc:
            logger.warning("過頂預報傳播失敗 cat=%s: %s", cat, exc)
            continue

        sl = obs["sl"]; cl = obs["cl"]; so = obs["so"]; co = obs["co"]
        x0 = obs["x0"]; y0 = obs["y0"]; z0 = obs["z0"]

        cat_passes: list[dict] = []

        for si, nid in enumerate(nids):
            r_ti = r_raw[si]   # (N_times, 3)
            e_ti = e_raw[si]   # (N_times,)
            ok   = (e_ti == 0)
            if not ok.any():
                continue

            # Vectorized ECI → ECEF → ENU → elevation
            xe = cg_all * r_ti[:, 0] + sg_all * r_ti[:, 1]
            ye = -sg_all * r_ti[:, 0] + cg_all * r_ti[:, 1]
            ze = r_ti[:, 2]
            dx = xe - x0; dy = ye - y0; dz = ze - z0
            U_enu = cl * co * dx + cl * so * dy + sl * dz
            E_enu = -so * dx + co * dy
            N_enu = -sl * co * dx - sl * so * dy + cl * dz
            rng   = np.sqrt(E_enu ** 2 + N_enu ** 2 + U_enu ** 2)
            safe  = np.where(rng > 0.001, rng, 0.001)
            el    = np.where(ok, np.rad2deg(np.arcsin(np.clip(U_enu / safe, -1.0, 1.0))), -90.0)

            above = el >= mask_deg
            if not above.any():
                continue

            transitions = np.diff(above.astype(int))
            rise_list   = list(np.where(transitions == 1)[0] + 1)
            set_list    = list(np.where(transitions == -1)[0] + 1)

            if above[0]:
                rise_list = [0] + rise_list
            if above[-1]:
                set_list  = set_list + [n_steps - 1]

            for ri, si_ in zip(rise_list, set_list):
                if ri >= si_:
                    continue
                seg_el     = el[ri: si_ + 1]
                mx_offset  = int(np.argmax(seg_el))
                max_el     = float(seg_el[mx_offset])
                mx_i       = ri + mx_offset
                t_rise     = times[ri]
                t_set      = times[min(si_, n_steps - 1)]
                t_max      = times[mx_i]
                duration_s = int((t_set - t_rise).total_seconds())
                cat_passes.append({
                    "norad_id":   nid,
                    "name":       idx[nid]["name"],
                    "t_rise_utc": t_rise.isoformat(),
                    "t_max_utc":  t_max.isoformat(),
                    "t_set_utc":  t_set.isoformat(),
                    "max_el_deg": round(max_el, 1),
                    "duration_s": duration_s,
                    "color":      cfg["color"],
                })

        cat_passes.sort(key=lambda x: x["t_rise_utc"])
        all_passes[cat] = cat_passes[:max_per_cat]

    result: dict[str, Any] = {}
    for cat, cfg in _OVERPASS_CATS.items():
        result[cat] = {
            "label":    cfg["label"],
            "sublabel": cfg["sublabel"],
            "color":    cfg["color"],
            "passes":   all_passes[cat],
        }

    return {
        "categories": result,
        "hours":      hours,
        "mask_deg":   mask_deg,
        "timestamp":  t0.isoformat(),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 時間軸：歷史/未來 TLE 索引 + 時間指定覆蓋/過頂計算
# ─────────────────────────────────────────────────────────────────────────────
_TIMELINE_DAYS = 30       # ±30 天


def _get_index_for_time(ts: datetime) -> dict[int, dict[str, Any]]:
    """
    回傳適用於 ts 時刻的衛星 TLE 索引。
    - ts >= now-1h  → 直接用現行索引（最快）
    - ts 在過去     → 查詢 DB 中 epoch_utc <= ts 的最新 TLE
    """
    now = datetime.now(timezone.utc)
    if ts >= now - timedelta(hours=1):
        return get_sat_index()

    db = _resolve_db()
    if db is None:
        return get_sat_index()

    # 確保 ts 為 UTC aware
    if ts.tzinfo is None:
        ts = ts.replace(tzinfo=timezone.utc)
    ts_str = ts.strftime("%Y-%m-%d %H:%M:%S")

    logger.info("歷史 TLE 查詢：%s", ts_str)
    try:
        with duckdb.connect(str(db), read_only=True) as con:
            sql = _tle_select_sql(con, extra_where=f"r.epoch_utc <= TIMESTAMP '{ts_str}'")
            rows = con.execute(sql).fetchall()
    except Exception as exc:
        logger.warning("歷史 TLE 查詢失敗 (%s)，退回現行索引: %s", ts_str, exc)
        return get_sat_index()

    if not rows:
        logger.warning("歷史 TLE 查詢無結果 (%s)", ts_str)
        return get_sat_index()

    csv_meta = load_sat_metadata_csv()
    idx: dict[int, dict[str, Any]] = {}
    for norad_id, raw_name, l1, l2, db_src, db_launch, db_intl in rows:
        nid  = int(norad_id)
        name = (raw_name or "").strip().lstrip("0 ") or f"OBJECT {nid}"
        ov   = csv_meta.get(nid, {})
        final_name = ov.get("name_en")     or name
        final_src  = ov.get("source_code") or db_src
        final_intl = ov.get("intl_code")   or db_intl
        csv_date_str = ov.get("launch_date", "")
        if csv_date_str:
            try:
                final_launch: datetime | None = datetime.strptime(
                    csv_date_str, "%Y-%m-%d").replace(tzinfo=timezone.utc)
            except ValueError:
                final_launch = db_launch
        else:
            final_launch = db_launch
        idx[nid] = {
            "name":          final_name,
            "line1":         l1.strip() if l1 else "",
            "line2":         l2.strip() if l2 else "",
            "country":       classify_country(final_src),
            "purpose":       ov.get("purpose") or classify_purpose(final_name),
            "era":           classify_era(final_launch, final_intl),
            "constellation": ov.get("constellation") or classify_constellation(final_name),
        }
    logger.info("歷史 TLE 索引完成: %d 筆 (%s)", len(idx), ts_str)
    return idx


def compute_taipei_coverage_at(ts: datetime, mask_deg: float = _MASK_DEG) -> dict[str, Any]:
    """Coverage at an arbitrary time ts (past or future)."""
    if ts.tzinfo is None:
        ts = ts.replace(tzinfo=timezone.utc)
    now = datetime.now(timezone.utc)

    idx = _get_index_for_time(ts)
    candidates = _get_overpass_candidates(idx)
    jd0, fr0 = jday(ts.year, ts.month, ts.day, ts.hour, ts.minute,
                    ts.second + ts.microsecond * 1e-6)
    obs = _TAIPEI_OBS

    all_nids = [nid for nids in candidates.values() for nid in nids]
    empty_result: dict[str, Any] = {
        "categories": {
            cat: {"label": cfg["label"], "sublabel": cfg["sublabel"],
                  "color": cfg["color"], "count": 0, "visible_count": 0, "satellites": []}
            for cat, cfg in _OVERPASS_CATS.items()
        },
        "timestamp":     ts.isoformat(),
        "mask_deg":      mask_deg,
        "is_historical": ts < now - timedelta(hours=1),
        "is_future":     ts > now + timedelta(hours=1),
    }
    if not all_nids:
        return empty_result

    line1s = [idx[n]["line1"] for n in all_nids]
    line2s = [idx[n]["line2"] for n in all_nids]
    err_arr, r_arr = _sgp4_propagate_raw(all_nids, line1s, line2s, ts)

    valid = err_arr == 0
    el_all = np.full(len(all_nids), -90.0)
    az_all = np.zeros(len(all_nids))
    rng_all = np.zeros(len(all_nids))
    llh_all = np.full((len(all_nids), 3), np.nan)

    if valid.any():
        el_v, az_v, rng_v = _eci_to_elaz(r_arr[valid], jd0, fr0, obs)
        el_all[valid]  = el_v
        az_all[valid]  = az_v
        rng_all[valid] = rng_v
        llh_all[valid] = eci_to_llh_batch(r_arr[valid], ts)

    nid_to_i = {nid: i for i, nid in enumerate(all_nids)}
    result_cats: dict[str, Any] = {}

    for cat, cfg in _OVERPASS_CATS.items():
        sats = []
        for nid in candidates[cat]:
            i = nid_to_i[nid]
            if err_arr[i] != 0 or np.isnan(llh_all[i, 0]):
                continue
            lat, lon, alt = llh_all[i, 0], llh_all[i, 1], llh_all[i, 2]
            if not (-500.0 < float(alt) < 80_000.0):
                continue
            el  = float(el_all[i])
            sats.append({
                "norad_id":  nid,
                "name":      idx[nid]["name"],
                "lat":       round(float(lat), 4),
                "lon":       round(float(lon), 4),
                "alt_km":    round(float(alt), 1),
                "el_deg":    round(el, 2),
                "az_deg":    round(float(az_all[i]), 2),
                "range_km":  round(float(rng_all[i]), 1),
                "visible":   el >= mask_deg,
                "color":     cfg["color"],
            })
        result_cats[cat] = {
            "label":         cfg["label"],
            "sublabel":      cfg["sublabel"],
            "color":         cfg["color"],
            "count":         len(sats),
            "visible_count": sum(1 for s in sats if s["visible"]),
            "satellites":    sats,
        }

    return {
        "categories":    result_cats,
        "timestamp":     ts.isoformat(),
        "mask_deg":      mask_deg,
        "is_historical": ts < now - timedelta(hours=1),
        "is_future":     ts > now + timedelta(hours=1),
    }


def predict_taipei_passes_at(
    ts:          datetime,
    hours:       float = 24.0,
    step_sec:    float = 60.0,
    mask_deg:    float = _MASK_DEG,
    max_per_cat: int   = 20,
) -> dict[str, Any]:
    """Pass predictions starting from an arbitrary ts."""
    if ts.tzinfo is None:
        ts = ts.replace(tzinfo=timezone.utc)
    now = datetime.now(timezone.utc)

    idx = _get_index_for_time(ts)
    candidates = _get_overpass_candidates(idx)
    obs = _TAIPEI_OBS

    n_steps = int(hours * 3600 / step_sec)
    times   = [ts + timedelta(seconds=i * step_sec) for i in range(n_steps)]
    jd_fr   = np.array([
        jday(tt.year, tt.month, tt.day, tt.hour, tt.minute,
             tt.second + tt.microsecond * 1e-6)
        for tt in times
    ])
    jds = np.ascontiguousarray(jd_fr[:, 0]); frs = np.ascontiguousarray(jd_fr[:, 1])

    T_cent   = ((jds - 2451545.0) + frs) / 36525.0
    gmst_all = np.deg2rad(
        (280.46061837 + 360.98564736629 * (jds - 2451545.0 + frs)
         + 0.000387933 * T_cent ** 2) % 360.0)
    cg_all = np.cos(gmst_all); sg_all = np.sin(gmst_all)

    all_passes: dict[str, list] = {cat: [] for cat in _OVERPASS_CATS}

    for cat, cfg in _OVERPASS_CATS.items():
        nids = candidates[cat]
        if not nids or not _HAS_SATREC_ARRAY:
            continue

        line1s = [idx[n]["line1"] for n in nids]
        line2s = [idx[n]["line2"] for n in nids]

        try:
            sa = _SatrecArray([Satrec.twoline2rv(l1, l2) for l1, l2 in zip(line1s, line2s)])
            e_raw, r_raw, _ = sa.sgp4(jds, frs)
        except Exception as exc:
            logger.warning("時間軸過頂預報失敗 cat=%s: %s", cat, exc)
            continue

        sl = obs["sl"]; cl = obs["cl"]; so = obs["so"]; co = obs["co"]
        x0 = obs["x0"]; y0 = obs["y0"]; z0 = obs["z0"]
        cat_passes: list[dict] = []

        for si, nid in enumerate(nids):
            r_ti = r_raw[si]; e_ti = e_raw[si]; ok = (e_ti == 0)
            if not ok.any():
                continue
            xe = cg_all * r_ti[:, 0] + sg_all * r_ti[:, 1]
            ye = -sg_all * r_ti[:, 0] + cg_all * r_ti[:, 1]
            ze = r_ti[:, 2]
            dx = xe - x0; dy = ye - y0; dz = ze - z0
            U  = cl * co * dx + cl * so * dy + sl * dz
            E  = -so * dx + co * dy
            N  = -sl * co * dx - sl * so * dy + cl * dz
            rng = np.sqrt(E**2 + N**2 + U**2)
            safe = np.where(rng > 0.001, rng, 0.001)
            el  = np.where(ok, np.rad2deg(np.arcsin(np.clip(U / safe, -1.0, 1.0))), -90.0)

            above = el >= mask_deg
            if not above.any():
                continue

            trans = np.diff(above.astype(int))
            rises = list(np.where(trans == 1)[0] + 1)
            sets  = list(np.where(trans == -1)[0] + 1)
            if above[0]:  rises = [0] + rises
            if above[-1]: sets  = sets + [n_steps - 1]

            for ri, si_ in zip(rises, sets):
                if ri >= si_:
                    continue
                seg = el[ri: si_ + 1]
                mx  = int(np.argmax(seg))
                cat_passes.append({
                    "norad_id":   nid,
                    "name":       idx[nid]["name"],
                    "t_rise_utc": times[ri].isoformat(),
                    "t_max_utc":  times[ri + mx].isoformat(),
                    "t_set_utc":  times[min(si_, n_steps - 1)].isoformat(),
                    "max_el_deg": round(float(seg[mx]), 1),
                    "duration_s": int((times[min(si_, n_steps - 1)] - times[ri]).total_seconds()),
                    "color":      cfg["color"],
                })

            if len(cat_passes) >= max_per_cat:
                break

        cat_passes.sort(key=lambda x: x["t_rise_utc"])
        all_passes[cat] = cat_passes[:max_per_cat]

    return {
        "categories": {
            cat: {"label": cfg["label"], "sublabel": cfg["sublabel"],
                  "color": cfg["color"], "passes": all_passes[cat]}
            for cat, cfg in _OVERPASS_CATS.items()
        },
        "hours":       hours,
        "mask_deg":    mask_deg,
        "timestamp":   ts.isoformat(),
        "is_historical": ts < now - timedelta(hours=1),
        "is_future":     ts > now + timedelta(hours=1),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 2-A  即時近距離配對掃描
# ─────────────────────────────────────────────────────────────────────────────
_conj_cache:     dict[str, Any] | None = None
_conj_loaded_at: float = 0.0


def build_conjunction_summary(
    threshold_km: float = _CONJ_THRESHOLD_KM,
    max_pairs:    int   = 200,
) -> dict[str, Any]:
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

    ok      = (err_arr == 0)
    ok_nids = [all_nids[i] for i in range(len(all_nids)) if ok[i]]
    ok_r    = r_arr[ok]

    if len(ok_nids) < 2:
        return {"count": 0, "pairs": [], "threshold_km": threshold_km,
                "total_scanned": len(ok_nids)}

    ok_llh = eci_to_llh_batch(ok_r, t)
    alt_ok = (ok_llh[:, 2] > -500.0) & (ok_llh[:, 2] < 80_000.0)
    filt_nids = [ok_nids[i] for i in range(len(ok_nids)) if alt_ok[i]]
    filt_r    = ok_r[alt_ok]
    filt_llh  = ok_llh[alt_ok]

    if len(filt_nids) < 2:
        return {"count": 0, "pairs": [], "threshold_km": threshold_km,
                "total_scanned": len(filt_nids)}

    t1        = time.monotonic()
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
# 嵌入式前端（主 3D 地球儀頁面）
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
#taipei-link{font-size:11px;color:#00E5FF;text-decoration:none;
  border:1px solid #00E5FF44;padding:3px 8px;border-radius:3px;white-space:nowrap;flex-shrink:0}
#taipei-link:hover{background:#00E5FF22}
#main{display:flex;flex:1;overflow:hidden}
#panel{width:300px;min-width:260px;background:#0d1117;
  border-right:1px solid #21262d;display:flex;flex-direction:column;overflow:hidden}
#search-box{display:flex;align-items:center;gap:4px;
  padding:6px 10px 4px;border-bottom:1px solid #21262d;flex-shrink:0}
#search-input{flex:1;background:#161b22;border:1px solid #21262d;border-radius:4px;
  color:#c9d1d9;font-size:11px;padding:5px 8px;outline:none;transition:.15s}
#search-input:focus{border-color:#58a6ff}
#search-input::placeholder{color:#484f58}
#search-clear{background:none;border:none;color:#484f58;cursor:pointer;
  font-size:15px;line-height:1;padding:0 2px;transition:.12s;display:none}
#search-clear:hover{color:#c9d1d9}
#tabs{display:flex;border-bottom:1px solid #21262d;flex-shrink:0}
.tab{flex:1;padding:6px 0;font-size:11px;color:#8b949e;text-align:center;
  cursor:pointer;border:none;background:transparent;transition:.15s;border-bottom:2px solid transparent}
.tab:hover{color:#c9d1d9;background:#161b22}
.tab.active{color:#58a6ff;border-bottom:2px solid #58a6ff}
#panel-back{display:none;padding:5px 10px 4px;border-bottom:1px solid #21262d;flex-shrink:0}
#panel-back button{background:none;border:none;color:#8b949e;font-size:11px;cursor:pointer;
  display:flex;align-items:center;gap:5px}
#panel-back button:hover{color:#c9d1d9}
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
.conj-row{padding:6px 10px;border-bottom:1px solid #161b22;cursor:pointer;transition:.12s}
.conj-row:hover{background:#1c2128}
.conj-km{font-size:11px;font-weight:bold;color:#FF6B6B;min-width:50px}
.conj-names{flex:1;font-size:11px;color:#c9d1d9;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.conj-note{font-size:9px;color:#484f58;padding:4px 10px 8px;font-style:italic}
#filter-status{padding:5px 10px;font-size:11px;color:#8b949e;
  border-top:1px solid #21262d;flex-shrink:0;min-height:26px}
#layer-bar{border-top:1px solid #21262d;background:#0d1117;padding:5px 0 3px;flex-shrink:0}
#layer-bar .sec{padding:3px 12px 2px;font-size:10px;color:#484f58;text-transform:uppercase;letter-spacing:.8px}
#layer-bar label{display:flex;align-items:center;gap:7px;padding:4px 12px;cursor:pointer;
  color:#8b949e;font-size:11px;transition:.12s;user-select:none}
#layer-bar label:hover{background:#1c2128;color:#c9d1d9}
#layer-bar input[type=checkbox]{margin:0;cursor:pointer;accent-color:#58a6ff}
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
#db-info-bar{display:flex;align-items:center;gap:12px;padding:2px 14px;
  font-size:10px;background:#050810;border-bottom:1px solid #161b22;
  flex-shrink:0;overflow-x:auto;white-space:nowrap;color:#8b949e}
.dbi-sep{color:#21262d;margin:0 2px}
.dbi-k{color:#484f58}
.dbi-v{color:#8b949e}
.dbi-v.ok{color:#3fb950}
.dbi-v.warn{color:#d29922}
</style>
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
  <a id="taipei-link" href="/taipei" target="_blank">&#127988; 台北覆蓋分析</a>
  <div id="ts">—</div>
</div>

<div id="db-info-bar">
  <span class="dbi-k">DB:</span>
  <span id="dbi-name" class="dbi-v">—</span>
  <span class="dbi-sep">|</span>
  <span class="dbi-k">更新：</span><span id="dbi-mtime" class="dbi-v">—</span>
  <span class="dbi-sep">|</span>
  <span class="dbi-k">有效衛星：</span><span id="dbi-sats" class="dbi-v">—</span>
  <span class="dbi-sep">|</span>
  <span class="dbi-k">TLE 日期：</span><span id="dbi-epoch" class="dbi-v">—</span>
  <span class="dbi-sep">|</span>
  <span class="dbi-k">大小：</span><span id="dbi-size" class="dbi-v">—</span>
</div>

<div id="main">
  <div id="panel">
    <div id="search-box">
      <input id="search-input" type="text" placeholder="搜尋 NORAD ID 或衛星名稱..." autocomplete="off"/>
      <button id="search-clear" onclick="clearSearch()">×</button>
    </div>
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

function setPanelMode(mode){
  panelMode=mode;
  const showTabs=(mode==='tabs');
  document.getElementById('tabs').style.display=showTabs?'flex':'none';
  document.getElementById('panel-back').style.display=showTabs?'none':'block';
  if(showTabs) renderPanel(activeTab);
}

function backToTabs(){
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

async function loadDbInfo(){
  try{
    const r=await fetch('/api/db_info');
    if(!r.ok) return;
    const d=await r.json();
    if(d.error) return;
    document.getElementById('dbi-name').textContent=d.db_name||'—';
    if(d.db_updated_at){
      const mt=new Date(d.db_updated_at);
      const diffH=(Date.now()-mt)/3600000;
      const el=document.getElementById('dbi-mtime');
      el.textContent=mt.toISOString().replace('T',' ').slice(0,16)+' UTC';
      el.className='dbi-v '+(diffH<48?'ok':'warn');
    }
    if(d.valid_sat_count!=null)
      document.getElementById('dbi-sats').textContent=d.valid_sat_count.toLocaleString()+' 顆';
    if(d.epoch_min&&d.epoch_max)
      document.getElementById('dbi-epoch').textContent=d.epoch_min.slice(0,10)+' ~ '+d.epoch_max.slice(0,10);
    if(d.db_size_mb!=null)
      document.getElementById('dbi-size').textContent=d.db_size_mb.toLocaleString()+' MB';
  }catch(e){console.warn('DB info 載入失敗',e);}
}

async function init(){
  await initCesium();
  await loadBordersLayer();
  await loadStats();
  loadConjunctions();
  loadDbInfo();
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
# 台北覆蓋分析 2D 頁面 HTML 模板（Token 由路由替換）
# ─────────────────────────────────────────────────────────────────────────────
_TAIPEI_HTML_TMPL = r"""<!DOCTYPE html>
<html lang="zh-TW">
<head>
<meta charset="utf-8"/>
<title>台北覆蓋分析 — Cesium 2D</title>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<style>
*{box-sizing:border-box;margin:0;padding:0}
html,body{height:100%;overflow:hidden;background:#0a0e17;color:#c9d1d9;font-family:Tahoma,sans-serif}
#topbar{display:flex;align-items:center;gap:10px;padding:6px 14px;
  background:linear-gradient(90deg,#0d1117,#161b22);border-bottom:1px solid #1e3a5f;
  height:44px;flex-shrink:0}
#topbar h1{font-size:13px;color:#58a6ff;white-space:nowrap;flex:1}
.tbtn{font-size:11px;background:#162032;border:1px solid #1e3a5f;color:#8b949e;
  padding:3px 10px;border-radius:3px;cursor:pointer;transition:.15s;text-decoration:none;
  display:inline-flex;align-items:center}
.tbtn:hover{color:#c9d1d9;border-color:#58a6ff}
#clock{font-size:11px;color:#8b949e;white-space:nowrap}
#main{display:flex;height:calc(100vh - 44px)}
#panel{width:290px;background:#0d1117;border-right:1px solid #21262d;
  display:flex;flex-direction:column;overflow:hidden;flex-shrink:0}
#panel-header{padding:8px 12px 5px;border-bottom:1px solid #21262d;flex-shrink:0}
#panel-header h2{font-size:12px;color:#58a6ff;margin-bottom:2px}
#panel-header .sub{font-size:10px;color:#484f58}
#tab-bar{display:flex;border-bottom:1px solid #21262d;flex-shrink:0}
.ptab{flex:1;padding:5px 0;font-size:11px;color:#8b949e;text-align:center;
  cursor:pointer;background:none;border:none;border-bottom:2px solid transparent;transition:.15s}
.ptab:hover{color:#c9d1d9;background:#161b22}
.ptab.active{color:#58a6ff;border-bottom-color:#58a6ff}
#panel-body{flex:1;overflow-y:auto;padding:4px 0}
#panel-body::-webkit-scrollbar{width:3px}
#panel-body::-webkit-scrollbar-thumb{background:#21262d}
.cat-card{margin:6px 8px;border-radius:5px;background:#111827;
  border:1px solid #21262d;overflow:hidden;cursor:pointer;transition:.15s}
.cat-card:hover{border-color:#30404d}
.cat-card.active{border-color:#58a6ff;background:#162032}
.cat-header{display:flex;align-items:center;gap:8px;padding:8px 10px 4px}
.cat-dot{width:10px;height:10px;border-radius:50%;flex-shrink:0}
.cat-label{flex:1;font-size:12px;font-weight:bold;color:#e6edf3}
.cat-cnt{font-size:11px;color:#8b949e}
.cat-sublabel{font-size:10px;color:#484f58;padding:0 10px 5px}
.cat-stats{display:flex;gap:5px;padding:0 8px 8px}
.cstat{flex:1;text-align:center;background:#0d1117;border-radius:3px;padding:4px 2px}
.cstat .sv{font-size:15px;font-weight:bold}
.cstat .sl{font-size:9px;color:#484f58;text-transform:uppercase;letter-spacing:.5px}
.pass-item{padding:7px 10px;border-bottom:1px solid #161b22}
.pass-name{font-size:11px;font-weight:bold;color:#c9d1d9;white-space:nowrap;
  overflow:hidden;text-overflow:ellipsis;margin-bottom:3px}
.pass-row{display:flex;align-items:center;gap:8px}
.pass-time{font-size:10px;color:#8b949e}
.pass-el{font-size:11px;font-weight:bold}
.pass-dur{font-size:9px;color:#484f58;margin-top:2px}
.pass-cat-hdr{padding:6px 10px 4px;background:#161b22;font-size:10px;
  font-weight:bold;text-transform:uppercase;letter-spacing:.6px;
  display:flex;align-items:center;gap:6px;border-bottom:1px solid #21262d}
.pass-empty{padding:10px 12px;font-size:11px;color:#484f58;font-style:italic}
#map-wrap{flex:1;position:relative}
#cesiumContainer{position:absolute;inset:0}
#map-loading{position:absolute;inset:0;background:rgba(10,14,23,.85);
  display:flex;align-items:center;justify-content:center;z-index:20;
  font-size:13px;color:#58a6ff;letter-spacing:1px}
#legend{position:absolute;bottom:28px;right:10px;z-index:15;
  background:rgba(13,17,23,.9);border:1px solid #21262d;border-radius:5px;
  padding:8px 12px;min-width:170px}
#legend h4{font-size:10px;color:#8b949e;margin-bottom:6px;
  text-transform:uppercase;letter-spacing:.6px}
.leg-row{display:flex;align-items:center;gap:7px;margin-bottom:4px;
  cursor:pointer;border-radius:3px;padding:1px 3px;transition:.12s}
.leg-row:hover{background:#161b22}
.leg-row.muted{opacity:.4}
.leg-dot{width:9px;height:9px;border-radius:50%;flex-shrink:0}
.leg-label{font-size:11px;color:#c9d1d9}
#statusbar{position:absolute;bottom:0;left:0;right:0;height:22px;z-index:15;
  background:rgba(10,14,23,.8);display:flex;align-items:center;gap:12px;padding:0 10px}
.sitem{font-size:10px;color:#8b949e}
.sitem span{color:#58a6ff}
.cesium-widget-credits{display:none!important}
#timeline{position:absolute;bottom:22px;left:0;right:0;z-index:15;
  background:rgba(10,14,23,.92);border-top:1px solid #21262d;padding:6px 12px 5px}
#tl-slider{width:100%;accent-color:#58a6ff;cursor:pointer;height:4px;margin-bottom:4px}
#tl-row{display:flex;align-items:center;gap:8px}
#tl-bounds{font-size:9px;color:#484f58;white-space:nowrap}
#tl-time{font-size:11px;color:#e6edf3;flex:1;text-align:center;font-weight:bold}
.tl-now-btn{font-size:10px;background:#162032;border:1px solid #1e3a5f;color:#58a6ff;
  padding:2px 8px;border-radius:3px;cursor:pointer;transition:.15s;white-space:nowrap}
.tl-now-btn:hover{background:#1e3a5f}
.tl-badge{font-size:10px;padding:2px 7px;border-radius:3px;font-weight:bold;white-space:nowrap}
.tl-badge.now{color:#4CAF50;background:rgba(76,175,80,.15)}
.tl-badge.hist{color:#FF9800;background:rgba(255,152,0,.15)}
.tl-badge.future{color:#58a6ff;background:rgba(88,166,255,.15)}
#tl-loading{font-size:10px;color:#484f58}
</style>
<script>
(function(){
  function _css(h){var l=document.createElement('link');l.rel='stylesheet';l.href=h;document.head.appendChild(l);}
  function _js(src,ok,fail){var s=document.createElement('script');s.src=src;s.onload=ok;s.onerror=fail||function(){};document.head.appendChild(s);}
  function _err(m){var el=document.getElementById('map-loading');if(el){el.textContent=m;el.style.color='#f44336';}}
  _css('/cesium/Widgets/widgets.css');
  _js('/cesium/Cesium.js',function(){startApp();},function(){_err('本機 Cesium 未安裝，請將 CesiumJS 解壓至 data/cesium/');});
}());
</script>
</head>
<body>
<div id="topbar">
  <h1>&#127988; 台北覆蓋分析 &#8212; 2000 km 範圍衛星過頂</h1>
  <div id="clock">&#8212;</div>
  <button class="tbtn" onclick="refreshAll()">&#8635; 更新</button>
  <a class="tbtn" href="/">&#8592; 返回主頁</a>
</div>
<div id="main">
  <div id="panel">
    <div id="panel-header">
      <h2>衛星分類覆蓋</h2>
      <div class="sub">台北 25.03&#176;N 121.57&#176;E &nbsp;&#183;&nbsp; 仰角遮蔽 5&#176;</div>
    </div>
    <div id="tab-bar">
      <button class="ptab active" data-tab="overview" onclick="switchPanelTab('overview')">分類概況</button>
      <button class="ptab" data-tab="passes" onclick="switchPanelTab('passes')">過頂預報</button>
    </div>
    <div id="panel-body"><div class="pass-empty">載入中...</div></div>
  </div>
  <div id="map-wrap">
    <div id="cesiumContainer"></div>
    <div id="map-loading">初始化地圖...</div>
    <div id="legend">
      <h4>衛星分類</h4>
      <div class="leg-row" data-cat="US_EO" onclick="toggleCatFilter('US_EO')">
        <div class="leg-dot" style="background:#4488FF"></div>
        <span class="leg-label">美國商用光學</span>
      </div>
      <div class="leg-row" data-cat="CN_COMM" onclick="toggleCatFilter('CN_COMM')">
        <div class="leg-dot" style="background:#FF9800"></div>
        <span class="leg-label">中國商用光學</span>
      </div>
      <div class="leg-row" data-cat="CN_MIL" onclick="toggleCatFilter('CN_MIL')">
        <div class="leg-dot" style="background:#F44336"></div>
        <span class="leg-label">中國軍用偵察</span>
      </div>
      <div class="leg-row" data-cat="TW_TASA" onclick="toggleCatFilter('TW_TASA')">
        <div class="leg-dot" style="background:#00E5FF"></div>
        <span class="leg-label">台灣 TASA</span>
      </div>
      <div style="border-top:1px solid #21262d;margin-top:5px;padding-top:5px;
                  display:flex;align-items:center;gap:7px">
        <div style="width:14px;height:3px;background:#FFD600;border-radius:1px;flex-shrink:0"></div>
        <span style="font-size:10px;color:#8b949e">2000 km 覆蓋圈</span>
      </div>
    </div>
    <!-- Timeline bar -->
    <div id="timeline">
      <input type="range" id="tl-slider" min="-43200" max="43200" step="15" value="0"
             oninput="tlSliderInput(this.value)" onchange="tlSliderChange(this.value)"/>
      <div id="tl-row">
        <span id="tl-bounds" class="tl-lb">&#8722;30天</span>
        <span id="tl-time">&#8212;</span>
        <span id="tl-bounds" class="tl-rb">&#43;30天</span>
        <button class="tl-now-btn" onclick="jumpToNow()">&#9654; 現在</button>
        <span id="tl-badge" class="tl-badge now">現在</span>
        <span id="tl-loading"></span>
      </div>
    </div>
    <div id="statusbar">
      <div class="sitem">衛星 <span id="st-total">&#8212;</span></div>
      <div class="sitem">可見 <span id="st-visible">&#8212;</span></div>
      <div class="sitem">時刻 <span id="st-time">&#8212;</span></div>
    </div>
  </div>
</div>
<script>
'use strict';
function startApp(){

const TAIPEI_LAT=25.0330, TAIPEI_LON=121.5654;
const CATS={
  US_EO:  {label:'美國商用光學衛星',sublabel:'Vantor/Maxar · Planet SkySat',color:'#4488FF'},
  CN_COMM:{label:'中國商用光學衛星',sublabel:'SuperView · 高分 · 吉林',  color:'#FF9800'},
  CN_MIL: {label:'中國軍用偵察衛星',sublabel:'遙感 Yaogan',                        color:'#F44336'},
  TW_TASA:{label:'台灣 TASA 衛星',  sublabel:'Formosat-5 / -7 / -8',              color:'#00E5FF'},
};

let viewer=null, satDs=null, circleDs=null;
let coverageData=null, passesData=null;
let activePanelTab='overview';
let activeCatFilter=null;

// ── Timeline state ─────────────────────────────────────────────────────────
let _tlMin=0;            // minutes offset from now (0=now, negative=past)
let _tlDebounce=null;
let _autoTimer=null;
let _coverageCtrl=null;   // AbortController for in-flight coverage fetch
let _passesCtrl=null;     // AbortController for in-flight passes fetch
let _loading=false;       // prevents concurrent _loadForTs calls
let _lastAction=0;        // rate-limit timestamp for buttons
const _TL_MAX=43200;     // 30 days in minutes

// ── Clock ──────────────────────────────────────────────────────────────────
function _cst(d){
  return new Date(d.getTime()+8*3600*1000).toISOString().replace('T',' ').slice(0,19)+' CST';
}
document.getElementById('clock').textContent=_cst(new Date());
setInterval(function(){document.getElementById('clock').textContent=_cst(new Date());},1000);

// ── Cesium init ────────────────────────────────────────────────────────────
async function initCesium(){
  Cesium.Ion.defaultAccessToken='CESIUM_TOKEN_PLACEHOLDER';

  viewer=new Cesium.Viewer('cesiumContainer',{
    // Initialise directly in 2D — no morph transition, no drift
    sceneMode:           Cesium.SceneMode.SCENE2D,
    mapProjection:       new Cesium.WebMercatorProjection(),
    animation:           false,
    timeline:            false,
    baseLayerPicker:     false,
    // OpenStreetMap tiles — good regional detail at zoom 5-8
    imageryProvider: new Cesium.UrlTemplateImageryProvider({
      url:'https://tile.openstreetmap.org/{z}/{x}/{y}.png',
      maximumLevel:18,
      credit:'© OpenStreetMap contributors',
    }),
    terrainProvider:     new Cesium.EllipsoidTerrainProvider(),
    sceneModePicker:     false,
    infoBox:             false,
    geocoder:            false,
    homeButton:          false,
    navigationHelpButton:false,
    selectionIndicator:  false,
    fullscreenButton:    false,
  });

  viewer.cesiumWidget.creditContainer.style.display='none';
  // Freeze clock — prevents time-driven camera updates that cause 2D drift
  viewer.clock.shouldAnimate=false;

  // setView() is synchronous — unlike flyTo() it does NOT start an animation
  // that can drift the camera position after load.
  viewer.camera.setView({
    destination:Cesium.Rectangle.fromDegrees(100.0,9.0,144.0,41.0),
  });

  // Disable rotation and tilt — meaningless in 2D and can cause drift
  viewer.scene.screenSpaceCameraController.enableRotate=false;
  viewer.scene.screenSpaceCameraController.enableTilt=false;

  satDs=new Cesium.CustomDataSource('satellites');
  circleDs=new Cesium.CustomDataSource('circle');
  await viewer.dataSources.add(satDs);
  await viewer.dataSources.add(circleDs);
  _drawCircle();

  document.getElementById('map-loading').style.display='none';
}

function _drawCircle(){
  const R=6371.0, cover=2000.0, ar=cover/R, n=128;
  const lat0=Cesium.Math.toRadians(TAIPEI_LAT);
  const lon0=Cesium.Math.toRadians(TAIPEI_LON);
  const pts=[];
  for(let i=0;i<=n;i++){
    const b=(i/n)*2*Math.PI;
    const la=Math.asin(Math.sin(lat0)*Math.cos(ar)+Math.cos(lat0)*Math.sin(ar)*Math.cos(b));
    const lo=lon0+Math.atan2(Math.sin(b)*Math.sin(ar)*Math.cos(lat0),
                             Math.cos(ar)-Math.sin(lat0)*Math.sin(la));
    pts.push(Cesium.Cartesian3.fromRadians(lo,la));
  }
  circleDs.entities.add({polyline:{positions:pts,width:2,
    material:new Cesium.ColorMaterialProperty(Cesium.Color.fromCssColorString('#FFD600').withAlpha(0.75)),
    arcType:Cesium.ArcType.NONE}});
  circleDs.entities.add({position:Cesium.Cartesian3.fromDegrees(TAIPEI_LON,TAIPEI_LAT),
    point:{pixelSize:9,color:Cesium.Color.fromCssColorString('#FFD600'),
           outlineColor:Cesium.Color.BLACK,outlineWidth:1.5}});
}

// ── Data loading ───────────────────────────────────────────────────────────
async function loadCoverage(ts=null){
  if(_coverageCtrl) _coverageCtrl.abort();
  _coverageCtrl=new AbortController();
  const url=ts?'/api/taipei_coverage_at?ts='+encodeURIComponent(ts):'/api/taipei_coverage';
  try{
    const r=await fetch(url,{signal:_coverageCtrl.signal});
    if(!r.ok) throw new Error('HTTP '+r.status);
    coverageData=await r.json();
    _updateMarkers();
    _updateStatus();
    if(activePanelTab==='overview') renderOverview();
  }catch(e){
    if(e.name==='AbortError') return;
    console.warn('Coverage error',e);
  }
}

async function loadPasses(ts=null){
  if(_passesCtrl) _passesCtrl.abort();
  _passesCtrl=new AbortController();
  const url=ts?'/api/taipei_passes_at?ts='+encodeURIComponent(ts):'/api/taipei_passes';
  try{
    const r=await fetch(url,{signal:_passesCtrl.signal});
    if(!r.ok) throw new Error('HTTP '+r.status);
    passesData=await r.json();
    if(activePanelTab==='passes') renderPasses();
  }catch(e){
    if(e.name==='AbortError') return;
    console.warn('Passes error',e);
  }
}

function refreshAll(){
  const now=Date.now();
  if(now-_lastAction<1000) return;
  _lastAction=now;
  clearTimeout(_tlDebounce);
  clearTimeout(_autoTimer);
  _loadForTs();
}

// ── Markers ────────────────────────────────────────────────────────────────
function _updateMarkers(){
  if(!coverageData||!satDs) return;
  satDs.entities.suspendEvents();
  satDs.entities.removeAll();
  Object.entries(coverageData.categories).forEach(([catId,cd])=>{
    if(activeCatFilter&&activeCatFilter!==catId) return;
    const col=Cesium.Color.fromCssColorString(cd.color);
    (cd.satellites||[]).forEach(s=>{
      satDs.entities.add({
        id:'sat_'+s.norad_id,
        position:Cesium.Cartesian3.fromDegrees(s.lon,s.lat,0),
        point:{
          pixelSize:s.visible?8:5,
          color:s.visible?col:col.withAlpha(0.4),
          outlineColor:s.visible?Cesium.Color.WHITE.withAlpha(0.85):Cesium.Color.WHITE.withAlpha(0.2),
          outlineWidth:s.visible?1.5:0.5,
        },
        label:{
          text:s.name,
          font:'11px Tahoma,sans-serif',
          fillColor:Cesium.Color.WHITE,
          outlineColor:Cesium.Color.BLACK,
          outlineWidth:2,
          style:Cesium.LabelStyle.FILL_AND_OUTLINE,
          pixelOffset:new Cesium.Cartesian2(0,-14),
          horizontalOrigin:Cesium.HorizontalOrigin.CENTER,
          verticalOrigin:Cesium.VerticalOrigin.BOTTOM,
          distanceDisplayCondition:new Cesium.DistanceDisplayCondition(0,4.6e6),
          disableDepthTestDistance:Number.POSITIVE_INFINITY,
        },
      });
    });
  });
  satDs.entities.resumeEvents();
}

// ── Category filter ────────────────────────────────────────────────────────
function toggleCatFilter(catId){
  activeCatFilter=(activeCatFilter===catId)?null:catId;
  document.querySelectorAll('.leg-row[data-cat]').forEach(el=>{
    el.classList.toggle('muted',activeCatFilter!==null&&el.dataset.cat!==activeCatFilter);
  });
  document.querySelectorAll('.cat-card').forEach(el=>{
    el.classList.toggle('active',activeCatFilter!==null&&el.dataset.cat===activeCatFilter);
  });
  _updateMarkers();
}

// ── Panel tabs ─────────────────────────────────────────────────────────────
function switchPanelTab(tab){
  activePanelTab=tab;
  document.querySelectorAll('.ptab').forEach(el=>{
    el.classList.toggle('active',el.dataset.tab===tab);
  });
  if(tab==='overview') renderOverview();
  else renderPasses();
}

// ── Overview panel ─────────────────────────────────────────────────────────
function renderOverview(){
  const body=document.getElementById('panel-body');
  body.innerHTML='';
  if(!coverageData){body.innerHTML="<div class='pass-empty'>資料載入中...</div>";return;}
  Object.entries(CATS).forEach(([catId,cfg])=>{
    const cd=coverageData.categories[catId];
    if(!cd) return;
    const card=document.createElement('div');
    card.className='cat-card'+(activeCatFilter===catId?' active':'');
    card.dataset.cat=catId;
    card.innerHTML=
      "<div class='cat-header'><div class='cat-dot' style='background:"+cfg.color+"'></div>"
      +"<span class='cat-label'>"+cfg.label+"</span>"
      +"<span class='cat-cnt'>"+cd.count+" 顆</span></div>"
      +"<div class='cat-sublabel'>"+cfg.sublabel+"</div>"
      +"<div class='cat-stats'>"
      +"<div class='cstat'><div class='sv' style='color:"+cfg.color+"'>"+cd.count+"</div><div class='sl'>資料庫</div></div>"
      +"<div class='cstat'><div class='sv' style='color:#4CAF50'>"+cd.visible_count+"</div><div class='sl'>可見</div></div>"
      +"</div>";
    card.addEventListener('click',()=>toggleCatFilter(catId));
    body.appendChild(card);
  });
}

// ── Pass panel ─────────────────────────────────────────────────────────────
function renderPasses(){
  const body=document.getElementById('panel-body');
  body.innerHTML='';
  if(!passesData){
    body.innerHTML="<div class='pass-empty'>過頂預報載入中（10–30 秒）...</div>";
    return;
  }
  const catIds=activeCatFilter?[activeCatFilter]:Object.keys(CATS);
  catIds.forEach(catId=>{
    const cfg=CATS[catId];
    const cd=passesData.categories[catId];
    if(!cd) return;
    const hdr=document.createElement('div');
    hdr.className='pass-cat-hdr';
    hdr.innerHTML="<span style='width:8px;height:8px;border-radius:50%;background:"+cfg.color+";display:inline-block'></span>"+cfg.label;
    body.appendChild(hdr);
    const passes=cd.passes||[];
    if(!passes.length){
      const e=document.createElement('div');
      e.className='pass-empty';
      e.textContent='24 小時內無過頂記錄';
      body.appendChild(e); return;
    }
    passes.slice(0,8).forEach(p=>{
      const item=document.createElement('div');
      item.className='pass-item';
      const rT=new Date(p.t_rise_utc);
      const elCol=p.max_el_deg>=45?'#4CAF50':p.max_el_deg>=20?'#FF9800':'#8b949e';
      const dm=Math.floor(p.duration_s/60), ds_=p.duration_s%60;
      const ts=new Date(rT.getTime()+8*3600*1000).toISOString().replace('T',' ').slice(5,16);
      item.innerHTML=
        "<div class='pass-name'>"+p.name+"</div>"
        +"<div class='pass-row'>"
        +"<span class='pass-time'>&#8599; "+ts+" CST</span>"
        +"<span class='pass-el' style='color:"+elCol+"'>Max "+p.max_el_deg+"&deg;</span>"
        +"</div>"
        +"<div class='pass-dur'>持續 "+dm+"m "+ds_+"s</div>";
      body.appendChild(item);
    });
  });
}

// ── Status bar ─────────────────────────────────────────────────────────────
function _updateStatus(){
  if(!coverageData) return;
  let tot=0,vis=0;
  Object.values(coverageData.categories).forEach(c=>{tot+=c.count||0;vis+=c.visible_count||0;});
  document.getElementById('st-total').textContent=tot;
  document.getElementById('st-visible').textContent=vis;
  const t=new Date(coverageData.timestamp);
  document.getElementById('st-time').textContent=
    new Date(t.getTime()+8*3600*1000).toISOString().replace('T',' ').slice(11,19)+' CST';
}

// ── Timeline ───────────────────────────────────────────────────────────────
function _tlTs(){
  // Returns ISO string for current slider offset, or null if "now"
  if(Math.abs(_tlMin)<5) return null;
  return new Date(Date.now()+_tlMin*60*1000).toISOString();
}

function _tlUpdateDisplay(){
  const ts=new Date(Date.now()+_tlMin*60*1000);
  document.getElementById('tl-time').textContent=_cst(ts);
  const badge=document.getElementById('tl-badge');
  if(Math.abs(_tlMin)<5){
    badge.textContent='現在'; badge.className='tl-badge now';
  } else if(_tlMin<0){
    badge.textContent='歷史'; badge.className='tl-badge hist';
  } else {
    badge.textContent='預測'; badge.className='tl-badge future';
  }
}

function tlSliderInput(v){
  // Live display update while dragging — no API call yet
  _tlMin=parseInt(v);
  _tlUpdateDisplay();
  clearTimeout(_tlDebounce);
}

function tlSliderChange(v){
  // Mouseup / touch end: trigger data load with 500ms debounce
  _tlMin=parseInt(v);
  _tlUpdateDisplay();
  clearTimeout(_tlDebounce);
  clearTimeout(_autoTimer);
  _tlDebounce=setTimeout(()=>_loadForTs(),500);
}

function jumpToNow(){
  const now=Date.now();
  if(now-_lastAction<1000) return;
  _lastAction=now;
  _tlMin=0;
  document.getElementById('tl-slider').value=0;
  _tlUpdateDisplay();
  clearTimeout(_tlDebounce);
  clearTimeout(_autoTimer);
  _loadForTs();
}

async function _loadForTs(){
  if(_loading) return;
  _loading=true;
  const ts=_tlTs();
  const ldEl=document.getElementById('tl-loading');
  if(ldEl) ldEl.textContent='載入中...';
  try{
    await Promise.all([loadCoverage(ts), loadPasses(ts)]);
  }finally{
    _loading=false;
  }
  if(ldEl) ldEl.textContent='';
  // Schedule next auto-refresh only when at "now"
  if(!ts){
    clearTimeout(_autoTimer);
    _autoTimer=setTimeout(()=>_loadForTs(),60000);
  }
}

// ── Borders layer (2D) ─────────────────────────────────────────────────────
async function loadBorders(){
  try{
    // GeoJsonDataSource handles 2D polygon rendering correctly when fill is
    // transparent — gives clean national border outlines without filled shapes.
    const ds=await Cesium.GeoJsonDataSource.load('/api/layers/borders',{
      stroke:      Cesium.Color.fromCssColorString('#FFD600').withAlpha(0.7),
      strokeWidth: 1.5,
      fill:        Cesium.Color.TRANSPARENT,
      markerSymbol:' ',
    });
    await viewer.dataSources.add(ds);
  }catch(e){
    console.warn('Borders load failed',e);
  }
}

// ── Init ───────────────────────────────────────────────────────────────────
async function init(){
  await initCesium();
  loadBorders();         // fire-and-forget
  _tlUpdateDisplay();    // init timeline display
  await loadCoverage();
  loadPasses();
  // Start 60-second auto-refresh for "now" mode
  _autoTimer=setTimeout(()=>_loadForTs(),60000);
}

window.tlSliderInput=tlSliderInput;
window.tlSliderChange=tlSliderChange;
window.jumpToNow=jumpToNow;

init().catch(e=>{
  document.getElementById('map-loading').textContent='初始化失敗: '+e.message;
  console.error(e);
});

window.switchPanelTab=switchPanelTab;
window.toggleCatFilter=toggleCatFilter;
window.refreshAll=refreshAll;

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

    # ── 台北覆蓋分析 2D 頁面 ──────────────────────────────────────────────────

    @app.get("/taipei")
    def taipei_page():
        html = _TAIPEI_HTML_TMPL.replace("CESIUM_TOKEN_PLACEHOLDER", CESIUM_ION_TOKEN)
        return html, 200, {"Content-Type": "text/html; charset=utf-8"}

    @app.get("/api/taipei_coverage")
    def api_taipei_coverage():
        try:
            mask = float(request.args.get("mask_deg", _MASK_DEG))
            mask = max(0.0, min(mask, 85.0))
        except ValueError:
            mask = _MASK_DEG
        t0 = time.monotonic()
        data = compute_taipei_coverage(mask_deg=mask)
        data["elapsed_sec"] = round(time.monotonic() - t0, 2)
        resp = make_response(_json.dumps(data, ensure_ascii=False).encode("utf-8"))
        resp.headers["Content-Type"]  = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = "public, max-age=60"
        return resp

    @app.get("/api/taipei_passes")
    def api_taipei_passes():
        try:
            hours = float(request.args.get("hours", 24.0))
            hours = max(1.0, min(hours, 72.0))
        except ValueError:
            hours = 24.0
        try:
            step = float(request.args.get("step_sec", 60.0))
            step = max(10.0, min(step, 300.0))
        except ValueError:
            step = 60.0
        try:
            mask = float(request.args.get("mask_deg", _MASK_DEG))
            mask = max(0.0, min(mask, 85.0))
        except ValueError:
            mask = _MASK_DEG
        t0 = time.monotonic()
        data = predict_taipei_passes(hours=hours, step_sec=step, mask_deg=mask)
        data["elapsed_sec"] = round(time.monotonic() - t0, 2)
        resp = make_response(_json.dumps(data, ensure_ascii=False).encode("utf-8"))
        resp.headers["Content-Type"]  = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = "public, max-age=300"
        return resp

    # ── 時間軸：指定時刻覆蓋 ────────────────────────────────────────────────

    @app.get("/api/taipei_coverage_at")
    def api_taipei_coverage_at():
        ts_str = request.args.get("ts", "").strip()
        try:
            ts = datetime.fromisoformat(ts_str.replace("Z", "+00:00"))
            if ts.tzinfo is None:
                ts = ts.replace(tzinfo=timezone.utc)
        except (ValueError, TypeError):
            ts = datetime.now(timezone.utc)
        now = datetime.now(timezone.utc)
        ts  = max(ts, now - timedelta(days=_TIMELINE_DAYS))
        ts  = min(ts, now + timedelta(days=_TIMELINE_DAYS))
        try:
            mask = float(request.args.get("mask_deg", _MASK_DEG))
            mask = max(0.0, min(mask, 85.0))
        except ValueError:
            mask = _MASK_DEG
        t0   = time.monotonic()
        data = compute_taipei_coverage_at(ts, mask_deg=mask)
        data["elapsed_sec"] = round(time.monotonic() - t0, 2)
        resp = make_response(_json.dumps(data, ensure_ascii=False).encode("utf-8"))
        resp.headers["Content-Type"]  = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = "public, max-age=60"
        return resp

    @app.get("/api/taipei_passes_at")
    def api_taipei_passes_at():
        ts_str = request.args.get("ts", "").strip()
        try:
            ts = datetime.fromisoformat(ts_str.replace("Z", "+00:00"))
            if ts.tzinfo is None:
                ts = ts.replace(tzinfo=timezone.utc)
        except (ValueError, TypeError):
            ts = datetime.now(timezone.utc)
        now = datetime.now(timezone.utc)
        ts  = max(ts, now - timedelta(days=_TIMELINE_DAYS))
        ts  = min(ts, now + timedelta(days=_TIMELINE_DAYS))
        try:
            hours = float(request.args.get("hours", 24.0))
            hours = max(1.0, min(hours, 72.0))
        except ValueError:
            hours = 24.0
        try:
            step = float(request.args.get("step_sec", 60.0))
            step = max(10.0, min(step, 300.0))
        except ValueError:
            step = 60.0
        try:
            mask = float(request.args.get("mask_deg", _MASK_DEG))
            mask = max(0.0, min(mask, 85.0))
        except ValueError:
            mask = _MASK_DEG
        t0   = time.monotonic()
        data = predict_taipei_passes_at(ts, hours=hours, step_sec=step, mask_deg=mask)
        data["elapsed_sec"] = round(time.monotonic() - t0, 2)
        resp = make_response(_json.dumps(data, ensure_ascii=False).encode("utf-8"))
        resp.headers["Content-Type"]  = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = "public, max-age=300"
        return resp

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

    @app.get("/api/db_info")
    def api_db_info():
        resp = make_response(
            _json.dumps(get_db_info(), ensure_ascii=False, default=str).encode("utf-8"))
        resp.headers["Content-Type"]  = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = f"public, max-age={_DB_INFO_TTL}"
        return resp

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

    # ── 單顆衛星查詢 ─────────────────────────────────────────────────────────

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

        if q.isdigit():
            nid = int(q)
            if nid in idx:
                matches.append({"norad_id": nid, **idx[nid], "score": 0})

        for nid, info in idx.items():
            if q_up in info["name"].upper():
                if not any(m["norad_id"] == nid for m in matches):
                    matches.append({"norad_id": nid, **info, "score": 1})
            if len(matches) >= 60:
                break

        matches.sort(key=lambda x: (x["score"], x["name"]))
        top = matches[:20]

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
        "Scenario 04-Cesium-advanced02 (向量化 SGP4 + 近距離掃描 + 搜尋 + 台北覆蓋 + 時間軸) 啟動 — "
        "http://%s:%d  台北覆蓋(時間軸): http://%s:%d/taipei",
        HOST, PORT, HOST, PORT,
    )
    logger.info(
        "SatrecArray=%s  KD-tree=%s  接近閾值=%.0f km  快取 TTL=%d s",
        _HAS_SATREC_ARRAY, _HAS_KDTREE, _CONJ_THRESHOLD_KM, _CONJ_TTL,
    )
    logger.info("預熱衛星索引…")
    get_sat_index()
    get_stats()
    app.run(host=HOST, port=PORT, debug=True)
