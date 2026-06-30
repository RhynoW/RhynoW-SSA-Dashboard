#!/usr/bin/env python3
"""
scenario05-6.py
===============
PRC 衛星當前位置一鍵顯示系統（3D CesiumJS ION + 衛星詳細資料）
★ Himawari-8/9 Band 16（13.3 μm FAR-IR）即時影像疊加圖層 ★
★ 衛星自動更新間隔：Radio Buttons 選擇 60 秒（預設）/ 10 秒 / 1 秒 ★
★ 可勾選向量圖層：全球國界 / SSN 觀測站 / VATSIM FIR 邊界 ★
★ 畫面右下角顯示 ATRDC LOGO（360×360 px）★
★ 地球初始視角以台灣為中心 ★
★ 衛星分類選取時背景顯示亮橘色 ★

以 scenario03-1.py 為範本，修改：
  - 衛星位置更新頻率：60 秒 → 10 秒
  - 新增向量圖層 toggle：全球國界、主要城市、主要軍事基地
  - GET /api/layers/borders        → Natural Earth 110m 國界 GeoJSON（CDN proxy + 快取）
  - GET /api/layers/cities         → 主要城市 GeoJSON（hardcoded 50+ 城市）
  - GET /api/layers/military_bases → 主要軍事基地 GeoJSON（hardcoded ~60 基地）

環境變數：
    CESIUM_ION_TOKEN  CesiumJS ION 存取金鑰（https://ion.cesium.com）
    DB_PATH           DuckDB 資料庫路徑（預設同目錄 space_db.duckdb）
    PORT              Flask 連接埠（預設 5007）
    TEST_SATS         >0 僅取前 N 顆（開發測試用）
    TLE_CACHE_TTL     TLE 快取有效秒數（預設 600）

用法：
    $env:CESIUM_ION_TOKEN="your_token_here"
    python scenario05.py
"""

from __future__ import annotations

import logging
import os
import re
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import duckdb
import io as _io

import numpy as np
import requests
from dotenv import load_dotenv
from flask import Flask, jsonify, make_response, request, send_from_directory
from flask_cors import CORS
from sgp4.api import Satrec, jday

try:
    import markdown as _md
    _MD_OK = True
except ImportError:
    _md = None       # type: ignore[assignment]
    _MD_OK = False

try:
    from PIL import Image as _PILImage
    from pyproj import Transformer as _Transformer
    from scipy.ndimage import map_coordinates as _map_coords
    _REPROJECT_OK = True
except ImportError as _re:
    _REPROJECT_OK = False
    _PILImage = _Transformer = _map_coords = None  # type: ignore[assignment]

# ─────────────────────────────────────────────────────────────────────────────
# 常數
# ─────────────────────────────────────────────────────────────────────────────
R_EARTH_KM = 6378.137
F_EARTH    = 1 / 298.257223563
E2         = F_EARTH * (2 - F_EARTH)

RAW_TABLE       = "raw_tle_archive"
_SCRIPT_DIR       = Path(__file__).resolve().parent
DEFAULT_DB        = str(_SCRIPT_DIR / "space_db.duckdb")
PRC_PROFILE_DIR   = _SCRIPT_DIR / "prc_maneuver" / "sat_profiles"
_TEXTURE_DIR        = _SCRIPT_DIR / "data" / "textures"
_JS_DIR             = _SCRIPT_DIR / "data" / "js"
_CESIUM_LOCAL_DIR   = _SCRIPT_DIR / "data" / "cesium"
_HIMAWARI_CACHE_DIR = _SCRIPT_DIR / "data" / "himawari"
_HIMAWARI_LOCAL_FILE = _HIMAWARI_CACHE_DIR / "himawari_latest.png"

load_dotenv(_SCRIPT_DIR / ".env")

logging.basicConfig(
    level=os.getenv("LOG_LEVEL", "INFO").upper(),
    format="%(asctime)s %(levelname)s %(name)s %(message)s",
)
logger = logging.getLogger("scenario05_6")

if not _MD_OK:
    logger.warning("markdown 套件未安裝，profile 將以純文字顯示（pip install markdown）")

# ─────────────────────────────────────────────────────────────────────────────
# Himawari-8/9 NICT API 設定
# ─────────────────────────────────────────────────────────────────────────────
# ── Himawari 座標處理參考 ──────────────────────────────────────────────────────
# 方法來源：https://smythdesign.com/blog/georeferencing-himawari/
#
# NICT 1d/550 全圓盤 PNG（地靜止投影）：
#   URL:  https://himawari.asia/img/FULL_24h/B16/1d/550/{YYYY}/{MM}/{DD}/{HHMMSS}_0_0.png
#   最新時間: https://himawari.asia/img/D531106/latest.json
#
# 投影參數（GEOS — Geostationary Earth Orbit）：
#   +proj=geos +h=35785863 +ellps=WGS84 +lon_0=140.7 +sweep=y
#   Sweep 方向：y（Himawari/Meteosat 標準；GOES 用 x）
#
# 全圓盤邊界（GEOS 座標空間）：
#   11,000 px × 1,000 m/px ÷ 2 = ±5,500,000 m（X 與 Y 方向相同）
#   對應 gdal_translate 參數：-a_ullr -5500000 5500000 5500000 -5500000
#
# 重投影方式（無 GDAL）：
#   pyproj.Transformer(EPSG:4326 → GEOS) 做反向映射
#   scipy.ndimage.map_coordinates 做雙線性內插
#   輸出：1440×720 等經緯度 WGS84 PNG（含 Alpha，圓盤外透明）
# ──────────────────────────────────────────────────────────────────────────────
HIMAWARI_LATEST_URL  = "https://himawari.asia/img/D531106/latest.json"
HIMAWARI_B16_TMPL    = (
    "https://himawari.asia/img/FULL_24h/B16/1d/550"
    "/{year}/{month:02d}/{day:02d}/{hour:02d}{minute:02d}{second:02d}_0_0.png"
)
HIMAWARI_REQUESTS_TIMEOUT = 15

# GEOS 投影參數（smythdesign 方法）
_GEOS_CRS    = "+proj=geos +h=35785863 +ellps=WGS84 +lon_0=140.7 +sweep=y +no_defs"
_DISK_EXTENT = 5_500_000.0   # GEOS 空間半徑（公尺）；11000 px × 1000 m/px ÷ 2
_REPROJ_W    = 1440           # 輸出等經緯度影像寬（px）
_REPROJ_H    = 720            # 輸出等經緯度影像高（px）

# 重投影快取：{timestamp_str: reprojected_png_bytes}
_hw_reproj_cache: dict[str, bytes] = {}

# ─────────────────────────────────────────────────────────────────────────────
# 向量圖層資料
# ─────────────────────────────────────────────────────────────────────────────
_borders_cache: bytes | None = None  # Natural Earth 110m 國界 GeoJSON
_vatsim_cache:  bytes | None = None  # VATSpy FIR/UIR Boundaries GeoJSON

_SSN_STATIONS_GEOJSON = {
    "type": "FeatureCollection",
    "features": [
        # ── GEODSS 光學站 ──
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-106.6599,33.8172]},"properties":{"name":"GEODSS Socorro","type":"光學/電光","location":"White Sands Missile Range, NM, USA","status":"active","notes":"深空目標光學追蹤"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-156.2578,20.7088]},"properties":{"name":"GEODSS Maui (AMOS)","type":"光學/電光","location":"Haleakalā, Hawaii, USA","status":"active","notes":"深空目標光學追蹤"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[72.4522,-7.4117]},"properties":{"name":"GEODSS Diego Garcia","type":"光學/電光","location":"British Indian Ocean Territory","status":"active","notes":"深空目標光學追蹤"}},
        # ── 其他光學/電光系統 ──
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-156.2600,20.7100]},"properties":{"name":"MAUI Space Surveillance (MSSS)","type":"光學/電光","location":"Haleakalā, Hawaii, USA","status":"active","notes":"先進光電系統，多光譜成像"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-156.2565,20.7076]},"properties":{"name":"AEOS Telescope","type":"光學/電光","location":"Haleakalā, Hawaii, USA","status":"active","notes":"先進電光感測器（3.67m 望遠鏡）"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-14.40,-7.97]},"properties":{"name":"Ascension Range Radar","type":"雷達","location":"Ascension Island","status":"active","notes":"南大西洋遙測/追蹤站"}},
        # ── 主要雷達系統 ──
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-86.21,30.57]},"properties":{"name":"AN/FPS-85","type":"雷達","location":"Eglin AFB, Florida, USA","status":"active","notes":"SSN最大功率相控陣雷達；第20太空監視中隊操作"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[167.750,8.730]},"properties":{"name":"Space Fence (AN/FPS-133)","type":"雷達","location":"Kwajalein Atoll, Marshall Islands","status":"active","notes":"S波段相控陣雷達，2020年起作戰；可追蹤10cm以下碎片"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[174.14,52.72]},"properties":{"name":"AN/FPS-108 Cobra Dane","type":"雷達","location":"Shemya Island, Alaska, USA","status":"active","notes":"相控陣雷達，兼飛彈預警與太空追蹤"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[31.13,70.37]},"properties":{"name":"GLOBUS II","type":"雷達","location":"Vardø, Norway","status":"active","notes":"X波段碟形雷達；挪威情報局（NIS）操作，數據納入SSN"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-71.518,42.261]},"properties":{"name":"HUSIR (Haystack)","type":"雷達","location":"Westford, Massachusetts, USA","status":"active","notes":"超寬頻雷達，MIT林肯實驗室"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-71.520,42.263]},"properties":{"name":"HAX (Haystack Auxiliary)","type":"雷達","location":"Westford, Massachusetts, USA","status":"active","notes":"Haystack輔助雷達"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-71.87,42.28]},"properties":{"name":"Millstone Hill Radar","type":"雷達","location":"North Grafton, Massachusetts, USA","status":"active","notes":"MIT林肯實驗室追蹤雷達"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[167.755,8.735]},"properties":{"name":"ALTAIR","type":"雷達","location":"Kwajalein Atoll, Marshall Islands","status":"active","notes":"Reagan Test Site (RTS) 深空雷達"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-96.37,47.12]},"properties":{"name":"PARCS (AN/FPQ-16)","type":"雷達","location":"Cavalier AFS, North Dakota, USA","status":"active","notes":"飛彈預警/太空監視相控陣雷達；第10太空預警中隊"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[167.745,8.725]},"properties":{"name":"ALCOR","type":"雷達","location":"Kwajalein Atoll, Marshall Islands","status":"active","notes":"C波段成像雷達，Reagan Test Site (RTS)"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[167.760,8.720]},"properties":{"name":"TRADEX","type":"雷達","location":"Kwajalein Atoll, Marshall Islands","status":"active","notes":"L/S波段追蹤雷達，Reagan Test Site (RTS)"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[167.745,8.740]},"properties":{"name":"MMW (Millimeter Wave)","type":"雷達","location":"Kwajalein Atoll, Marshall Islands","status":"active","notes":"毫米波精細追蹤雷達，Reagan Test Site (RTS)"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[167.758,8.728]},"properties":{"name":"GBR-P","type":"雷達","location":"Kwajalein Atoll, Marshall Islands","status":"active","notes":"X波段地基雷達原型"}},
        # ── 飛彈預警/SSN協作感測器（BMEWS / PAVE PAWS）──
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-70.54,41.75]},"properties":{"name":"Cape Cod SFS (PAVE PAWS)","type":"飛彈預警/協作","location":"Bourne, Massachusetts, USA","status":"active","notes":"AN/FPS-123，UHF相控陣；第6太空預警中隊"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-121.35,39.14]},"properties":{"name":"Beale AFB (PAVE PAWS)","type":"飛彈預警/協作","location":"California, USA","status":"active","notes":"AN/FPS-123，UHF相控陣；第7太空預警中隊"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-149.19,64.30]},"properties":{"name":"Clear SFS (BMEWS)","type":"飛彈預警/協作","location":"Alaska, USA","status":"active","notes":"AN/FPS-120，UHF相控陣；第13太空預警中隊"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-0.67,54.36]},"properties":{"name":"RAF Fylingdales (BMEWS)","type":"飛彈預警/協作","location":"England, UK","status":"active","notes":"AN/FPS-126，三面相控陣，360°覆蓋；英國皇家空軍操作"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-68.78,76.53]},"properties":{"name":"Thule/Pituffik SB (BMEWS)","type":"飛彈預警/協作","location":"Greenland","status":"active","notes":"AN/FPS-120；兼為SCN衛星追控站"}},
        # ── 衛星追控網（SCN）站點 ──
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-71.43,42.97]},"properties":{"name":"New Boston AFS","type":"衛星追控","location":"New Hampshire, USA","status":"active","notes":"23rd Space Operations Squadron SCN站點"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-120.57,34.75]},"properties":{"name":"Vandenberg SFB (SCN)","type":"衛星追控","location":"California, USA","status":"active","notes":"23rd Space Operations Squadron SCN站點"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-158.13,21.36]},"properties":{"name":"Kaena Point","type":"衛星追控","location":"Oahu, Hawaii, USA","status":"active","notes":"23rd Space Operations Squadron SCN站點"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[144.87,13.58]},"properties":{"name":"Guam Remote Tracking Station","type":"衛星追控","location":"Guam","status":"active","notes":"23rd Space Operations Squadron SCN站點"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-1.08,51.07]},"properties":{"name":"RAF Oakhanger","type":"衛星追控","location":"Hampshire, UK","status":"active","notes":"23rd Space Operations Squadron SCN站點（英國）"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[72.4508,-7.413]},"properties":{"name":"Diego Garcia (SCN)","type":"衛星追控","location":"British Indian Ocean Territory","status":"active","notes":"23rd Space Operations Squadron SCN站點（兼GEODSS）"}},
        # ── 數據中心 ──
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-104.80,38.70]},"properties":{"name":"Space Surveillance Center (SSC)","type":"數據中心","location":"Cheyenne Mountain Complex, Colorado, USA","status":"active","notes":"SPACECOM太空監視作戰中心"}},
        # ── 已除役站點 ──
        {"type":"Feature","geometry":{"type":"Point","coordinates":[40.04,37.95]},"properties":{"name":"AN/FPS-79 Pirinclik","type":"已除役","location":"Diyarbakır, Turkey","status":"decommissioned","notes":"原SSN輔助雷達，已關閉"}},
        {"type":"Feature","geometry":{"type":"Point","coordinates":[-5.609,37.17]},"properties":{"name":"MOSS (Morón)","type":"已除役","location":"Morón Air Base, Spain","status":"decommissioned","notes":"光學站，1997–2012年運作"}},
    ],
}


def _reproject_himawari_to_wgs84(img_bytes: bytes) -> bytes:
    """
    將 Himawari GEOS 全圓盤 PNG 重投影至 WGS84 等經緯度（-180~180, -90~90）。

    算法（反向映射）：
      1. 建立輸出格網：(lon, lat) ∈ [-180,180] × [-90,90]，1440×720 px
      2. pyproj EPSG:4326 → GEOS 反向轉換：(lon,lat) → (x_geos, y_geos)
      3. GEOS 座標轉原圖像素：px = (x/extent + 0.5) × src_w
      4. scipy map_coordinates 雙線性取樣
      5. 圓盤外（|x|>extent 或 inf）設 alpha=0（透明）
    """
    if not _REPROJECT_OK:
        return img_bytes  # fallback: 回傳原圖

    # ── 載入原圖 ──
    src_img = _PILImage.open(_io.BytesIO(img_bytes)).convert("RGBA")
    src_arr = np.array(src_img, dtype=np.float32)  # (H, W, 4)
    src_h, src_w = src_arr.shape[:2]

    # ── 輸出格網：等經緯度，Y 軸由北至南 ──
    lons = np.linspace(-180.0, 180.0, _REPROJ_W, endpoint=False) + 180.0 / _REPROJ_W
    lats = np.linspace( 90.0, -90.0, _REPROJ_H, endpoint=False) -  90.0 / _REPROJ_H
    lon_grid, lat_grid = np.meshgrid(lons, lats)   # (720, 1440)

    # ── EPSG:4326 → GEOS（反向映射：目標格網 → 來源投影空間）──
    trans = _Transformer.from_crs("EPSG:4326", _GEOS_CRS, always_xy=True)
    flat_lons = lon_grid.ravel()
    flat_lats = lat_grid.ravel()
    x_geos, y_geos = trans.transform(flat_lons, flat_lats, errcheck=False)
    x_geos = np.array(x_geos, dtype=np.float64).reshape(_REPROJ_H, _REPROJ_W)
    y_geos = np.array(y_geos, dtype=np.float64).reshape(_REPROJ_H, _REPROJ_W)

    # ── 可見遮罩：在圓盤內且無 inf/nan ──
    valid = (
        np.isfinite(x_geos) & np.isfinite(y_geos)
        & (np.abs(x_geos) <= _DISK_EXTENT)
        & (np.abs(y_geos) <= _DISK_EXTENT)
    )

    # ── GEOS 座標 → 來源像素座標 ──
    # X 正方向向右，Y 正方向向上（影像 Y=0 在頂端，故取反）
    src_px_x = ( x_geos / (2.0 * _DISK_EXTENT) + 0.5) * src_w
    src_px_y = (-y_geos / (2.0 * _DISK_EXTENT) + 0.5) * src_h
    np.clip(src_px_x, 0, src_w - 1, out=src_px_x)
    np.clip(src_px_y, 0, src_h - 1, out=src_px_y)

    coords = np.array([src_px_y.ravel(), src_px_x.ravel()])  # (row, col)

    # ── 雙線性取樣（每通道）──
    out_arr = np.zeros((_REPROJ_H, _REPROJ_W, 4), dtype=np.uint8)
    for ch in range(4):
        sampled = _map_coords(src_arr[:, :, ch], coords, order=1, mode="nearest")
        out_arr[:, :, ch] = sampled.reshape(_REPROJ_H, _REPROJ_W).astype(np.uint8)

    # ── 圓盤外透明 ──
    out_arr[~valid, 3] = 0

    # ── 輸出 PNG ──
    out_img = _PILImage.fromarray(out_arr, mode="RGBA")
    buf = _io.BytesIO()
    out_img.save(buf, format="PNG", optimize=False)
    return buf.getvalue()


def _fetch_himawari_latest() -> dict[str, Any]:
    """
    向 NICT 查詢最新 Himawari 時間戳記。
    回傳 dict: {timestamp: str, year, month, day, hour, minute, second}
    若失敗則回傳 None。
    """
    try:
        resp = requests.get(HIMAWARI_LATEST_URL, timeout=HIMAWARI_REQUESTS_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
        # 典型回應: {"date": "2024-01-15 02:00:00", "file": "..."}
        date_str: str = data.get("date", "")
        if not date_str:
            logger.warning("Himawari latest.json 無 date 欄位: %s", data)
            return {}
        # 格式: "YYYY-MM-DD HH:MM:SS"
        dt = datetime.strptime(date_str.strip(), "%Y-%m-%d %H:%M:%S")
        return {
            "timestamp": date_str.strip(),
            "year":   dt.year,
            "month":  dt.month,
            "day":    dt.day,
            "hour":   dt.hour,
            "minute": dt.minute,
            "second": dt.second,
        }
    except Exception as exc:
        logger.warning("Himawari latest.json 取得失敗: %s", exc)
        return {}


def _build_b16_url(ts: dict[str, Any]) -> str:
    """根據時間戳記 dict 組出 Band 16 圖像 URL。"""
    return HIMAWARI_B16_TMPL.format(**ts)


# ─────────────────────────────────────────────────────────────────────────────
# Cesium infoBox 樣式
# ─────────────────────────────────────────────────────────────────────────────
_INFOBOX_CSS = (
    "<style>"
    "body{font-family:Tahoma,sans-serif;font-size:13px;color:#e8e8e8;background:#1e1e1e;margin:0;padding:2px 8px 10px}"
    "h1{font-size:14px;color:#8ab4f8;margin:6px 0 8px;padding-bottom:4px;border-bottom:2px solid #3a5a8a}"
    "h2{font-size:12px;color:#90caf9;margin:12px 0 4px;padding-bottom:2px;"
    "   border-bottom:1px solid #2a3a4a;text-transform:uppercase;letter-spacing:.5px}"
    "table{border-collapse:collapse;width:100%;margin:4px 0;font-size:12px}"
    "td,th{padding:3px 8px;border:1px solid #3a3a3a;vertical-align:top;color:#e8e8e8}"
    "th{background:#2a3a4a;font-weight:600;text-align:left;white-space:nowrap}"
    "tr:nth-child(even) td{background:#252525}"
    "p{margin:5px 0;line-height:1.55}"
    "blockquote{background:#2a2a1e;border-left:3px solid #ffc107;"
    "           margin:6px 0;padding:4px 10px;border-radius:0 4px 4px 0;font-size:12px}"
    "ul,ol{margin:4px 0;padding-left:20px}"
    "li{margin:2px 0;font-size:12px}"
    "a{color:#8ab4f8;text-decoration:none}"
    "a:hover{text-decoration:underline}"
    "hr{border:none;border-top:1px solid #3a3a3a;margin:8px 0}"
    "strong{font-weight:600;color:#e8e8e8}"
    "</style>"
)

# ─────────────────────────────────────────────────────────────────────────────
# 衛星分類規則
# ─────────────────────────────────────────────────────────────────────────────
CATEGORY_RULES: list[tuple[str, list[str]]] = [
    ("商業光學",      ["JILIN", "SUPERVIEW"]),
    ("軍事/偵察",     ["YAOGAN", "GJZ"]),
    ("軍事光學/SAR",  ["GAOFEN", "TIANHUI"]),
    ("氣象/環境",     ["FENGYUN", "YUNHAI", "DAQI", "HJ-", "HJ "]),
    ("商業通訊/IoT",  ["GEESAT", "TIANQI", "CENTISPACE", "QIANFAN", "DEAR",
                       "SHIKONGXING", "XIWANG", "HONGHU", "HJS"]),
    ("空間站補給",    ["TIANZHOU"]),
    ("試驗衛星",      ["SHIYAN", "SY-", "SJ-", "TIANZHI", "PRC TEST",
                       "HXMT", "HUIYAN", "PRC_TEST"]),
    ("火箭殘骸",      ["CZ-2C R/B", "CZ-2C RB", " R/B"]),
    ("不明物體",      ["OBJECT"]),
]

CATEGORY_COLORS: dict[str, str] = {
    "商業光學":      "#FFC107",
    "軍事/偵察":     "#F44336",
    "軍事光學/SAR":  "#FF5722",
    "氣象/環境":     "#03A9F4",
    "商業通訊/IoT":  "#4CAF50",
    "空間站補給":    "#E040FB",
    "試驗衛星":      "#9E9E9E",
    "火箭殘骸":      "#607D8B",
    "不明物體":      "#BDBDBD",
    "其他":          "#78909C",
}


def infer_category(name: str) -> str:
    n = name.upper()
    for cat, keywords in CATEGORY_RULES:
        if any(k.upper() in n for k in keywords):
            return cat
    return "其他"


def load_prc_catalog() -> dict[int, dict[str, Any]]:
    """讀取 sat_profiles/*.md，解析 NORAD ID、衛星名稱、分類。"""
    catalog: dict[int, dict[str, Any]] = {}
    if not PRC_PROFILE_DIR.exists():
        logger.warning("sat_profiles 目錄不存在: %s", PRC_PROFILE_DIR)
        return catalog
    for md_file in PRC_PROFILE_DIR.glob("*.md"):
        m = re.match(r"^(\d+)_(.+)$", md_file.stem)
        if not m:
            continue
        norad_id = int(m.group(1))
        name = m.group(2).replace("_", " ")
        catalog[norad_id] = {
            "norad_id":    norad_id,
            "name":        name,
            "category":    infer_category(name),
            "profile_path": md_file,
        }
    logger.info("PRC 衛星清單載入完成: %d 顆", len(catalog))
    return catalog


# ─────────────────────────────────────────────────────────────────────────────
# 環境設定
# ─────────────────────────────────────────────────────────────────────────────
DB_PATH          = Path(os.getenv("DB_PATH", DEFAULT_DB))
HOST             = os.getenv("HOST", "0.0.0.0")
PORT             = int(os.getenv("PORT", "5007"))
TEST_SATS        = int(os.getenv("TEST_SATS", "0"))
TLE_CACHE_TTL    = int(os.getenv("TLE_CACHE_TTL", "600"))
CESIUM_ION_TOKEN = os.getenv("CESIUM_ION_TOKEN", "")

if not CESIUM_ION_TOKEN:
    logger.warning("CESIUM_ION_TOKEN 未設定，使用離線底圖")

_full_catalog = load_prc_catalog()

if TEST_SATS > 0:
    _kept = sorted(_full_catalog.keys())[:TEST_SATS]
    prc_catalog = {nid: _full_catalog[nid] for nid in _kept}
    logger.warning("【測試模式】衛星數限制為 %d 顆（全目錄共 %d 顆）",
                   TEST_SATS, len(_full_catalog))
else:
    prc_catalog = _full_catalog

_tle_cache: dict[int, dict[str, str]] = {}
_tle_cache_loaded_at: float = 0.0


def get_tle_cache() -> dict[int, dict[str, str]]:
    global _tle_cache, _tle_cache_loaded_at
    age = time.monotonic() - _tle_cache_loaded_at
    if not _tle_cache or age > TLE_CACHE_TTL:
        logger.info("重新載入 TLE 快取（%d 顆）…", len(prc_catalog))
        t0 = time.monotonic()
        _tle_cache = load_tle_batch(list(prc_catalog.keys()))
        _tle_cache_loaded_at = time.monotonic()
        logger.info("TLE 快取更新完成：取得 %d / %d 顆，耗時 %.2f 秒",
                    len(_tle_cache), len(prc_catalog), _tle_cache_loaded_at - t0)
    return _tle_cache


# ─────────────────────────────────────────────────────────────────────────────
# 軌道計算
# ─────────────────────────────────────────────────────────────────────────────
def eci_to_llh(r_eci: np.ndarray, t_utc: datetime) -> tuple[float, float, float]:
    x, y, z = r_eci
    jd, fr = jday(
        t_utc.year, t_utc.month, t_utc.day,
        t_utc.hour, t_utc.minute,
        t_utc.second + t_utc.microsecond * 1e-6,
    )
    T        = ((jd - 2451545.0) + fr) / 36525.0
    gmst     = 280.46061837 + 360.98564736629 * (jd - 2451545.0 + fr) + 0.000387933 * T ** 2
    gmst_rad = np.deg2rad(gmst % 360.0)
    x_ecef   = np.cos(gmst_rad) * x + np.sin(gmst_rad) * y
    y_ecef   = -np.sin(gmst_rad) * x + np.cos(gmst_rad) * y
    z_ecef   = z
    lon = np.arctan2(y_ecef, x_ecef)
    r   = np.sqrt(x_ecef ** 2 + y_ecef ** 2)
    lat = np.arctan2(z_ecef, r * (1 - E2))
    alt = 0.0
    for _ in range(5):
        sin_lat = np.sin(lat)
        N       = R_EARTH_KM / np.sqrt(1 - E2 * sin_lat ** 2)
        cos_lat = np.cos(lat)
        alt     = (r / cos_lat - N) if abs(cos_lat) > 1e-9 else (abs(z_ecef) / (1.0 - E2) - N)
        lat     = np.arctan2(z_ecef, r * (1 - E2 * (N / (N + alt))))
    return float(np.rad2deg(lat)), float(np.rad2deg(lon)), float(alt)


def propagate_to_now(line1: str, line2: str) -> dict[str, float] | None:
    try:
        sat    = Satrec.twoline2rv(line1, line2)
        now    = datetime.now(timezone.utc)
        jd, fr = jday(now.year, now.month, now.day,
                      now.hour, now.minute,
                      now.second + now.microsecond * 1e-6)
        err, r_eci, _ = sat.sgp4(jd, fr)
        if err != 0:
            return None
        lat, lon, alt_km = eci_to_llh(np.array(r_eci), now)
        if not (-500.0 < alt_km < 80_000.0):
            logger.debug("SGP4 高度異常 %.1f km，略過", alt_km)
            return None
        return {"lat": lat, "lon": lon, "alt_km": round(alt_km, 1)}
    except Exception as exc:
        logger.debug("SGP4 傳播錯誤: %s", exc)
        return None


def propagate_arc(line1: str, line2: str,
                  hours: float = 2.0, pts: int = 120) -> list[dict[str, float]]:
    sat    = Satrec.twoline2rv(line1, line2)
    now    = datetime.now(timezone.utc)
    step_s = hours * 3600.0 / pts
    positions: list[dict[str, float]] = []
    for i in range(pts + 1):
        t = now + timedelta(seconds=i * step_s)
        jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute,
                      t.second + t.microsecond * 1e-6)
        err, r_eci, _ = sat.sgp4(jd, fr)
        if err != 0:
            continue
        try:
            lat, lon, alt_km = eci_to_llh(np.array(r_eci), t)
        except Exception:
            continue
        if not (-500.0 < alt_km < 80_000.0):
            continue
        positions.append({"lat": round(lat, 4), "lon": round(lon, 4), "alt_km": round(alt_km, 1)})
    return positions


def _resolve_db() -> Path | None:
    if DB_PATH.exists():
        return DB_PATH
    alt = DB_PATH.parent / "space_db_slim.duckdb"
    if alt.exists():
        logger.info("使用 slim DB: %s", alt)
        return alt
    logger.error("找不到資料庫: %s 或 space_db_slim.duckdb", DB_PATH)
    return None


def load_tle_batch(norad_ids: list[int]) -> dict[int, dict[str, str]]:
    if not norad_ids:
        return {}
    db = _resolve_db()
    if db is None:
        return {}
    id_str = ",".join(map(str, norad_ids))
    try:
        with duckdb.connect(str(db), read_only=True) as con:
            rows = con.execute(f"""
                SELECT norad_id, line1, line2
                FROM {RAW_TABLE}
                WHERE norad_id IN ({id_str})
                  AND line1 IS NOT NULL
                  AND line2 IS NOT NULL
                QUALIFY ROW_NUMBER() OVER (
                    PARTITION BY norad_id ORDER BY epoch_utc DESC
                ) = 1
            """).fetchall()
    except Exception as exc:
        logger.error("資料庫查詢失敗: %s", exc)
        return {}
    return {int(r[0]): {"line1": r[1], "line2": r[2]} for r in rows}


# ─────────────────────────────────────────────────────────────────────────────
# 嵌入式前端頁面（CesiumJS 1.114 + Himawari B16 疊加層）
# ─────────────────────────────────────────────────────────────────────────────
_HTML = """<!DOCTYPE html>
<html lang="zh-TW">
<head>
<meta charset="utf-8"/>
<title>PRC 衛星軌道預測 — Scenario 05-6 (Himawari B16 + 向量圖層 + LOGO)</title>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<link href="https://cesium.com/downloads/cesiumjs/releases/1.114/Build/Cesium/Widgets/widgets.css" rel="stylesheet"/>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:'Segoe UI',Tahoma,sans-serif;background:#0d1117;color:#c9d1d9;height:100vh;display:flex;overflow:hidden}
#panel{width:240px;min-width:200px;background:#161b22;border-right:1px solid #21262d;display:flex;flex-direction:column;overflow:hidden;z-index:1}
#ph{padding:12px 14px 10px;background:#0d1117;border-bottom:1px solid #21262d}
#ph h1{font-size:13px;color:#58a6ff;margin-bottom:2px}
#ph .sub{font-size:10px;color:#6e7681}
#stat{padding:6px 14px;font-size:11px;color:#8b949e;background:#161b22;border-bottom:1px solid #21262d;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.row-btns{display:flex;gap:5px;padding:7px 10px;border-bottom:1px solid #21262d}
.mb{flex:1;padding:5px 0;border-radius:4px;border:1px solid #30363d;background:transparent;color:#8b949e;cursor:pointer;font-size:11px;transition:.15s}
.mb:hover{background:#21262d;color:#c9d1d9}
.sec{padding:6px 14px 3px;font-size:10px;color:#484f58;text-transform:uppercase;letter-spacing:.8px}
#cat-list{flex:1;overflow-y:auto}
#cat-list::-webkit-scrollbar{width:4px}
#cat-list::-webkit-scrollbar-track{background:#161b22}
#cat-list::-webkit-scrollbar-thumb{background:#30363d;border-radius:2px}
.cb{display:flex;align-items:center;gap:7px;padding:7px 14px;cursor:pointer;border:none;background:transparent;color:#8b949e;font-size:12px;width:100%;text-align:left;transition:.15s}
.cb:hover{background:#1c2128;color:#c9d1d9}
.cb.on{background:#e65100;color:#fff}
.cb.on:hover{background:#d84800}
.dot{width:9px;height:9px;border-radius:50%;flex-shrink:0}
.cnt{margin-left:auto;font-size:10px;color:#484f58;background:#21262d;padding:1px 5px;border-radius:6px}
.cb.on .cnt{color:#fff;background:rgba(0,0,0,0.28)}
#footer{padding:10px 14px;border-top:1px solid #21262d;font-size:10px;color:#6e7681}
#footer label{display:flex;align-items:center;gap:6px;cursor:pointer;margin-top:6px}
#ts{color:#8b949e}
/* ── Himawari 狀態列 ── */
#hw-bar{padding:7px 14px;border-top:1px solid #21262d;font-size:10px;color:#6e7681;background:#161b22}
#hw-bar label{display:flex;align-items:center;gap:6px;cursor:pointer;font-size:11px;color:#8b949e}
#hw-bar label:hover{color:#c9d1d9}
#hw-ts{font-size:10px;color:#484f58;margin-top:3px;padding-left:22px;min-height:12px}
#hw-status{font-size:10px;padding-left:22px;min-height:12px}
/* ── 底圖選擇列 ── */
#basemap-bar{border-top:1px solid #21262d;background:#161b22;padding-bottom:4px}
#basemap-bar .sec{padding:7px 14px 3px;font-size:10px;color:#484f58;text-transform:uppercase;letter-spacing:.8px}
.bm-btn{display:block;width:calc(100% - 28px);margin:3px 14px;padding:5px 10px;background:transparent;border:1px solid #30363d;border-radius:4px;color:#8b949e;font-size:11px;text-align:left;cursor:pointer;transition:.12s}
.bm-btn:hover{background:#1c2128;color:#c9d1d9;border-color:#58a6ff}
.bm-btn.active{background:#1c2128;color:#58a6ff;border-color:#58a6ff}
.bm-btn:disabled{opacity:.4;cursor:not-allowed}
/* ── 向量圖層列 ── */
#layer-bar{border-top:1px solid #21262d;background:#161b22;padding-bottom:4px}
#layer-bar .sec{padding:7px 14px 3px;font-size:10px;color:#484f58;text-transform:uppercase;letter-spacing:.8px}
#layer-bar .cb{display:flex;align-items:center;gap:7px;padding:6px 14px;cursor:pointer;color:#8b949e;font-size:12px;transition:.12s;user-select:none}
#layer-bar .cb:hover{background:#1c2128;color:#c9d1d9}
#layer-bar .cb input{margin:0;cursor:pointer;accent-color:#58a6ff}
#cesiumContainer{flex:1;position:relative}
#warn{position:absolute;top:10px;left:50%;transform:translateX(-50%);
  background:rgba(244,67,54,.88);color:#fff;padding:6px 16px;border-radius:4px;
  font-size:12px;z-index:10;pointer-events:none;display:none;white-space:nowrap}
.cesium-widget-credits{display:none!important}
.cesium-infoBox-description{background:#1e1e1e!important;color:#e8e8e8!important}
.cesium-infoBox-description a{color:#8ab4f8!important}
.cesium-infoBox-description table,.cesium-infoBox-description td,.cesium-infoBox-description th{color:#e8e8e8!important;border-color:#3a3a3a!important}
#logo-overlay{position:absolute;bottom:0;right:0;width:180px;height:180px;pointer-events:none;z-index:5}
#logo-overlay img{width:100%;height:100%;object-fit:contain;display:block}
#ar-radios{display:flex;flex-direction:column;gap:3px;margin-top:4px;padding:0 14px}
#ar-radios label{display:flex;align-items:center;gap:7px;cursor:pointer;font-size:11px;color:#8b949e;padding:2px 0}
#ar-radios label:hover{color:#c9d1d9}
#ar-radios input[type=radio]{margin:0;cursor:pointer;accent-color:#58a6ff}
</style>
</head>
<body>

<div id="panel">
  <div id="ph">
    <h1>&#127760; PRC 衛星軌道預測</h1>
    <div class="sub">Scenario 05 &middot; 3D CesiumJS &middot; Himawari B16 &middot; 向量圖層</div>
  </div>
  <div id="stat">正在初始化...</div>
  <div class="row-btns">
    <button class="mb" onclick="selectAll()">全選</button>
    <button class="mb" onclick="clearAll()">清除</button>
    <button class="mb" onclick="doRefresh()">更新</button>
  </div>
  <div class="sec">衛星分類</div>
  <div id="cat-list"></div>

  <!-- Himawari B16 控制列 -->
  <div id="hw-bar">
    <label>
      <input type="checkbox" id="hw-toggle" onchange="toggleHimawari(this)"/>
      <span>&#127777; Himawari B16 (13.3&#181;m)</span>
    </label>
    <div id="hw-ts"></div>
    <div id="hw-status"></div>
  </div>

  <!-- 底圖選擇 -->
  <div id="basemap-bar">
    <div class="sec">地球底圖</div>
  </div>

  <!-- 向量圖層 -->
  <div id="layer-bar">
    <div class="sec">向量圖層</div>
    <label class="cb">
      <input type="checkbox" id="chk-borders" onchange="toggleLayer('borders',this)"/>
      <span class="dot" style="background:#FFD600;width:9px;height:9px;border-radius:50%;flex-shrink:0"></span>
      <span>全球國界</span>
    </label>
    <label class="cb">
      <input type="checkbox" id="chk-ssn" onchange="toggleLayer('ssn',this)"/>
      <span class="dot" style="background:#00E5FF;width:9px;height:9px;border-radius:50%;flex-shrink:0"></span>
      <span>SSN 地面觀測站</span>
    </label>
    <label class="cb">
      <input type="checkbox" id="chk-vatsim" onchange="toggleLayer('vatsim',this)"/>
      <span class="dot" style="background:#80DEEA;width:9px;height:9px;border-radius:50%;flex-shrink:0"></span>
      <span>VATSIM FIR 邊界</span>
    </label>
  </div>

  <div id="footer">
    <div>&#128339; <span id="ts">&#8212;</span> UTC</div>
    <div style="margin:7px 14px 3px;font-size:10px;color:#484f58;text-transform:uppercase;letter-spacing:.8px">自動更新間隔</div>
    <div id="ar-radios">
      <label><input type="radio" name="ar-interval" value="60000" checked onchange="setUpdateInterval(60000)"/> 60 秒（預設）</label>
      <label><input type="radio" name="ar-interval" value="10000" onchange="setUpdateInterval(10000)"/> 10 秒</label>
      <label><input type="radio" name="ar-interval" value="1000"  onchange="setUpdateInterval(1000)"/> 1 秒</label>
    </div>
    <label style="margin-top:6px"><input type="checkbox" id="show-arc" checked/> 點擊顯示軌道弧（2 h）</label>
  </div>
</div>

<div id="cesiumContainer">
  <div id="warn">&#9888; CESIUM_ION_TOKEN 未設定，使用離線底圖</div>
  <div id="logo-overlay"><img src="/api/logo" alt="ATRDC"/></div>
</div>

<script src="https://cesium.com/downloads/cesiumjs/releases/1.114/Build/Cesium/Cesium.js"></script>
<script>
'use strict';

let cats=[], sel=new Set(), viewer=null, satDs=null, orbitEnt=null, arTimer=null;
const entMap=new Map();

// ── 向量圖層狀態 ──────────────────────────────────────────────────────────
let borderDs=null, ssnDs=null, vatsimDs=null;

// ── 底圖狀態 ──────────────────────────────────────────────────────────────
let _currentBasemap='default';
let _basemapMeta={};

// ── Himawari 圖層狀態 ─────────────────────────────────────────────────────
let hwLayer=null;        // Cesium ImageryLayer
let hwRefreshTimer=null; // setInterval handle
const HW_REFRESH_MS = 10 * 60 * 1000; // 10 分鐘

// Himawari 重投影後為標準 WGS84 等經緯度，覆蓋全球範圍
// 圓盤外像素已設為透明（alpha=0），Cesium 不會顯示非可視區域
const HW_BOUNDS = new Cesium.Rectangle(
  Cesium.Math.toRadians(-180.0),  // west
  Cesium.Math.toRadians(-90.0),   // south
  Cesium.Math.toRadians( 180.0),  // east
  Cesium.Math.toRadians( 90.0)   // north
);

// ── 底圖函式 ──────────────────────────────────────────────────────────────
function _buildProvider(key){
  const m=_basemapMeta[key]||{};
  const t=m.type||'single';
  if(t==='tms'){
    return new Cesium.TileMapServiceImageryProvider({
      url:m.tms_url||'',
      credit:m.credit||'',
    });
  }
  return new Cesium.SingleTileImageryProvider({
    url:`/api/globe_texture/${key}`,
    rectangle:Cesium.Rectangle.fromDegrees(-180,-90,180,90),
    credit:m.credit||'',
  });
}

function switchBasemap(key){
  if(key===_currentBasemap) return;
  _currentBasemap=key;
  localStorage.setItem('basemap',key);
  document.querySelectorAll('.bm-btn').forEach(b=>{
    b.classList.toggle('active',b.dataset.bk===key);
  });
  viewer.imageryLayers.removeAll();
  viewer.imageryLayers.addImageryProvider(_buildProvider(key));
}

async function initBasemapBar(){
  let list=[];
  try{
    const r=await fetch('/api/textures');
    if(r.ok) list=await r.json();
  }catch(e){ console.warn('底圖列表載入失敗',e); }

  const bar=document.getElementById('basemap-bar');
  list.forEach(item=>{
    _basemapMeta[item.key]=item;
    const btn=document.createElement('button');
    btn.className='bm-btn';
    btn.dataset.bk=item.key;
    btn.textContent=item.label;
    if(!item.available) btn.disabled=true;
    btn.onclick=()=>switchBasemap(item.key);
    bar.appendChild(btn);
  });

  // 從 localStorage 還原底圖
  const saved=localStorage.getItem('basemap');
  const target=saved&&_basemapMeta[saved]&&_basemapMeta[saved].available
    ? saved : 'default';
  switchBasemap(target);
}

// ── 初始化 ─────────────────────────────────────────────────────────────────
async function init(){
  let cfg={cesium_token:'',has_token:false};
  try{
    const r=await fetch('/api/config');
    if(r.ok) cfg=await r.json();
  }catch(e){ console.warn('設定載入失敗',e); }

  Cesium.Ion.defaultAccessToken = cfg.cesium_token || '';

  if(!cfg.has_token) document.getElementById('warn').style.display='block';

  const commonOpts={
    animation:false, timeline:false,
    baseLayerPicker:false,
    imageryProvider:new Cesium.SingleTileImageryProvider({
      url:'/api/globe_texture/default',
      rectangle:Cesium.Rectangle.fromDegrees(-180,-90,180,90),
      credit:'NASA Blue Marble',
    }),
    sceneModePicker:true, infoBox:true,
    geocoder:false, homeButton:true, navigationHelpButton:false,
    selectionIndicator:true,
  };

  let viewerOpts;
  if(cfg.has_token){
    viewerOpts={
      ...commonOpts,
      terrain:Cesium.Terrain.fromWorldTerrain({requestWaterMask:true}),
    };
  }else{
    viewerOpts={
      ...commonOpts,
      terrainProvider:new Cesium.EllipsoidTerrainProvider(),
    };
  }

  viewer=new Cesium.Viewer('cesiumContainer', viewerOpts);
  viewer.cesiumWidget.creditContainer.style.display='none';
  viewer.scene.globe.enableLighting=true;

  await initBasemapBar();

  // 初始視角：以台灣為中心（lon 121°E, lat 23.5°N），高度 5,000 km
  viewer.camera.setView({
    destination: Cesium.Cartesian3.fromDegrees(121.0, 23.5, 5000000),
    orientation: {
      heading: Cesium.Math.toRadians(0),
      pitch:   Cesium.Math.toRadians(-90),
      roll:    0.0,
    },
  });

  satDs=new Cesium.CustomDataSource('prc_satellites');
  await viewer.dataSources.add(satDs);

  viewer.selectedEntityChanged.addEventListener(ent=>{
    _clearOrbit();
    if(!ent||!ent.properties) return;
    const nid=ent.properties.norad_id?.getValue();
    if(!nid) return;
    if(document.getElementById('show-arc').checked){
      showOrbitArc(nid, ent.properties.color_hex?.getValue()||'#ffffff');
    }
    loadProfile(ent, nid);
  });

  await loadCats();
}

// ── Himawari B16 圖層管理 ─────────────────────────────────────────────────

function _hwSetStatus(msg, isErr){
  const el=document.getElementById('hw-status');
  el.textContent=msg;
  el.style.color=isErr?'#f44336':'#4caf50';
}

async function _buildHimawariLayer(){
  // 先查詢 latest timestamp
  let tsLabel='';
  try{
    const r=await fetch('/api/himawari/latest');
    if(r.ok){
      const d=await r.json();
      tsLabel=d.timestamp||'';
    }
  }catch(e){ /* 靜默失敗，直接用 frame 端點 */ }

  // 使用帶 cache-buster 的 frame 端點，避免瀏覽器快取舊圖
  const bust=Date.now();
  const provider=new Cesium.SingleTileImageryProvider({
    url:'/api/himawari/frame?_t='+bust,
    rectangle:HW_BOUNDS,
    credit:'Himawari-8/9 Band 16 (13.3μm) © NICT / JMA',
  });

  const layer=viewer.imageryLayers.addImageryProvider(provider);
  layer.alpha=0.65;
  layer.brightness=1.0;
  layer.contrast=1.2;

  if(tsLabel){
    document.getElementById('hw-ts').textContent='影像時間: '+tsLabel+' UTC';
  }
  _hwSetStatus('B16 圖層已載入', false);
  return layer;
}

async function toggleHimawari(cb){
  if(cb.checked){
    _hwSetStatus('載入中...', false);
    document.getElementById('hw-ts').textContent='';
    try{
      hwLayer=await _buildHimawariLayer();
      // 每 10 分鐘自動重整
      clearInterval(hwRefreshTimer);
      hwRefreshTimer=setInterval(refreshHimawari, HW_REFRESH_MS);
    }catch(e){
      _hwSetStatus('Himawari 資料無法取得', true);
      document.getElementById('hw-ts').textContent='';
      cb.checked=false;
      console.warn('Himawari 載入失敗', e);
    }
  }else{
    _removeHimawariLayer();
    _hwSetStatus('', false);
    document.getElementById('hw-ts').textContent='';
  }
}

function _removeHimawariLayer(){
  clearInterval(hwRefreshTimer);
  hwRefreshTimer=null;
  if(hwLayer){
    viewer.imageryLayers.remove(hwLayer, true);
    hwLayer=null;
  }
}

async function refreshHimawari(){
  if(!hwLayer) return;
  _hwSetStatus('重新整理中...', false);
  try{
    // 移除舊圖層，建立新圖層（SingleTileImageryProvider 無法熱更新 URL）
    viewer.imageryLayers.remove(hwLayer, true);
    hwLayer=null;
    hwLayer=await _buildHimawariLayer();
  }catch(e){
    _hwSetStatus('重新整理失敗', true);
    console.warn('Himawari 重新整理失敗', e);
  }
}

// ── Profile 載入 ───────────────────────────────────────────────────────────
async function loadProfile(ent, nid){
  ent.description=new Cesium.ConstantProperty(
    '<div style="padding:12px;font-family:sans-serif;color:#666;font-size:13px">&#9203; 載入衛星資料中&#8230;</div>'
  );
  try{
    const r=await fetch('/api/prc/profile/'+nid);
    if(!r.ok) throw new Error('HTTP '+r.status);
    const d=await r.json();
    if(viewer.selectedEntity===ent){
      ent.description=new Cesium.ConstantProperty(d.html);
    }
  }catch(e){
    if(viewer.selectedEntity===ent){
      ent.description=new Cesium.ConstantProperty(
        '<div style="padding:8px;font-family:sans-serif;color:#c62828;font-size:12px">&#9888; 資料載入失敗: '+e.message+'</div>'
      );
    }
    console.warn('Profile 載入失敗', e);
  }
}

// ── 分類面板 ───────────────────────────────────────────────────────────────
async function loadCats(){
  try{
    const r=await fetch('/api/prc/categories');
    if(!r.ok) throw new Error('HTTP '+r.status);
    cats=await r.json();
    if(!Array.isArray(cats)) throw new Error('分類回應格式錯誤');
    renderCats();
    document.getElementById('stat').textContent='請選擇分類以顯示衛星';
  }catch(e){
    document.getElementById('stat').textContent='分類載入失敗: '+e.message;
    console.error(e);
  }
}

function renderCats(){
  const el=document.getElementById('cat-list');
  el.innerHTML='';
  cats.forEach(c=>{
    const on=sel.has(c.category);
    const b=document.createElement('button');
    b.className='cb'+(on?' on':'');
    b.dataset.cat=c.category;
    b.innerHTML=
      '<span class="dot" style="background:'+c.color+'"></span>'+
      '<span>'+c.category+'</span>'+
      '<span class="cnt">'+c.count+'</span>';
    b.onclick=()=>toggle(c.category,b);
    el.appendChild(b);
  });
}

function toggle(cat,btn){
  sel.has(cat)?sel.delete(cat):sel.add(cat);
  btn.classList.toggle('on');
  fetchPos();
}

function selectAll(){ cats.forEach(c=>sel.add(c.category)); renderCats(); fetchPos(); }
function clearAll(){
  sel.clear();
  satDs.entities.removeAll(); entMap.clear();
  _clearOrbit();
  renderCats();
  document.getElementById('stat').textContent='未選擇任何分類';
}

// ── 位置更新 ───────────────────────────────────────────────────────────────
async function fetchPos(){
  if(!sel.size){ satDs.entities.removeAll(); entMap.clear(); return; }
  document.getElementById('stat').textContent='載入中...';
  try{
    const q=encodeURIComponent([...sel].join(','));
    const r=await fetch('/api/prc/positions?cats='+q);
    if(!r.ok) throw new Error('HTTP '+r.status);
    const d=await r.json();
    if(!d||!Array.isArray(d.satellites)) throw new Error('回應格式錯誤');
    renderEntities(d);
    const ts=new Date(d.timestamp);
    document.getElementById('ts').textContent=ts.toISOString().replace('T',' ').slice(0,19);
    const total=[...sel].reduce((s,c)=>{
      const x=cats.find(v=>v.category===c); return s+(x?x.count:0);
    },0);
    document.getElementById('stat').textContent='顯示 '+d.count+' 顆 / 共 '+total+' 顆';
  }catch(e){
    document.getElementById('stat').textContent='錯誤: '+e.message;
    console.error(e);
  }
}

function renderEntities(d){
  const incoming=new Set(d.satellites.map(s=>s.norad_id));
  entMap.forEach((ent,id)=>{
    if(!incoming.has(id)){ satDs.entities.remove(ent); entMap.delete(id); }
  });
  d.satellites.forEach(s=>{
    const pos=Cesium.Cartesian3.fromDegrees(s.lon, s.lat, s.alt_km*1000);
    const col=Cesium.Color.fromCssColorString(s.color);
    if(entMap.has(s.norad_id)){
      entMap.get(s.norad_id).position=new Cesium.ConstantPositionProperty(pos);
    }else{
      const ent=satDs.entities.add({
        id:'sat_'+s.norad_id,
        name:s.name,
        position:pos,
        point:{
          pixelSize:7,
          color:col,
          outlineColor:Cesium.Color.WHITE.withAlpha(0.35),
          outlineWidth:1,
          scaleByDistance:new Cesium.NearFarScalar(5e5,1.6,8e6,0.8),
        },
        label:{
          text:s.name,
          font:'11px "Segoe UI",sans-serif',
          fillColor:Cesium.Color.WHITE,
          outlineColor:Cesium.Color.BLACK,
          outlineWidth:2,
          style:Cesium.LabelStyle.FILL_AND_OUTLINE,
          pixelOffset:new Cesium.Cartesian2(0,-16),
          distanceDisplayCondition:new Cesium.DistanceDisplayCondition(0,2e6),
        },
        description:new Cesium.ConstantProperty(_buildBasicDesc(s)),
        properties:{
          norad_id:new Cesium.ConstantProperty(s.norad_id),
          color_hex:new Cesium.ConstantProperty(s.color),
        },
      });
      entMap.set(s.norad_id, ent);
    }
  });
}

function _buildBasicDesc(s){
  return '<div style="font-family:Tahoma,sans-serif;font-size:13px;padding:4px">'
    +'<p style="color:#1565c0;font-weight:bold;margin-bottom:6px">'+s.name+'</p>'
    +'<table style="border-collapse:collapse;color:#333">'
    +'<tr><td style="padding:2px 10px 2px 0;color:#555">NORAD</td><td><b>'+s.norad_id+'</b></td></tr>'
    +'<tr><td style="padding:2px 10px 2px 0;color:#555">分類</td><td>'+s.category+'</td></tr>'
    +'<tr><td style="padding:2px 10px 2px 0;color:#555">高度</td><td><b>'+s.alt_km+' km</b></td></tr>'
    +'<tr><td style="padding:2px 10px 2px 0;color:#555">緯度</td><td>'+s.lat.toFixed(3)+'&deg;</td></tr>'
    +'<tr><td style="padding:2px 10px 2px 0;color:#555">經度</td><td>'+s.lon.toFixed(3)+'&deg;</td></tr>'
    +'</table>'
    +'<p style="color:#9e9e9e;font-size:11px;margin-top:8px">點擊以載入完整衛星資料</p>'
    +'</div>';
}

// ── 軌道弧 ─────────────────────────────────────────────────────────────────
function _clearOrbit(){
  if(orbitEnt){ viewer.entities.remove(orbitEnt); orbitEnt=null; }
}

async function showOrbitArc(nid, colorHex){
  try{
    const r=await fetch('/api/prc/orbit?norad_id='+nid+'&hours=2&pts=120');
    if(!r.ok) return;
    const d=await r.json();
    if(!d.positions||d.positions.length<2) return;
    const coords=d.positions.flatMap(p=>[p.lon, p.lat, p.alt_km*1000]);
    const glowColor=Cesium.Color.fromCssColorString(colorHex||'#ffffff').withAlpha(0.9);
    orbitEnt=viewer.entities.add({
      polyline:{
        positions:Cesium.Cartesian3.fromDegreesArrayHeights(coords),
        width:1.8,
        material:new Cesium.PolylineGlowMaterialProperty({glowPower:0.1, color:glowColor}),
        clampToGround:false,
        arcType:Cesium.ArcType.NONE,
      },
    });
  }catch(e){ console.warn('軌道弧載入失敗',e); }
}

// ── SSN 站點類型顏色映射 ─────────────────────────────────────────────────
const SSN_TYPE_COLORS={
  '光學/電光':    '#FFC107',
  '雷達':         '#00E5FF',
  '飛彈預警/協作':'#FF7043',
  '衛星追控':     '#66BB6A',
  '數據中心':     '#CE93D8',
  '已除役':       '#78909C',
};

// ── 向量圖層管理 ─────────────────────────────────────────────────────────
async function toggleLayer(type, cb){
  const ref={
    borders:{get:()=>borderDs, set:v=>{borderDs=v;}, url:'/api/layers/borders'},
    ssn:    {get:()=>ssnDs,    set:v=>{ssnDs=v;},    url:'/api/layers/ssn_stations'},
    vatsim: {get:()=>vatsimDs, set:v=>{vatsimDs=v;}, url:'/api/layers/vatsim_boundaries'},
  }[type];
  if(!ref) return;

  if(cb.checked){
    try{
      let ds;
      if(type==='borders'){
        // ★ describe:()=>'' 在 GeoJsonDataSource 解析期間攔截所有 describe(properties)
        //   呼叫，阻止 Cesium enumerate properties →
        //   RangeError: Too many properties to enumerate
        ds=await Cesium.GeoJsonDataSource.load(ref.url,{
          describe:function(){ return ''; },
          clampToGround:false,
        });
        ds.entities.values.forEach(ent=>{
          if(ent.polygon){
            // arcType GEODESIC：大地測量弧線，避免奇怪的幾何渲染失真
            ent.polygon.arcType=new Cesium.ConstantProperty(Cesium.ArcType.GEODESIC);
            ent.polygon.fill=new Cesium.ConstantProperty(false);
            ent.polygon.outline=new Cesium.ConstantProperty(true);
            ent.polygon.outlineColor=new Cesium.ConstantProperty(
              Cesium.Color.fromCssColorString('#FFD600').withAlpha(0.75));
            ent.polygon.outlineWidth=new Cesium.ConstantProperty(1.5);
          }
          ent.description=undefined;
          if(ent.label) ent.label.show=new Cesium.ConstantProperty(false);
        });
      }else if(type==='ssn'){
        // SSN 站點：依 type 屬性著色；同樣傳入空 describe 防止 enumerate
        ds=await Cesium.GeoJsonDataSource.load(ref.url,{describe:function(){ return ''; }});
        // 先展開成靜態陣列，避免 forEach 迭代中呼叫 ds.entities.add() 干擾原始迭代器
        [...ds.entities.values].forEach(ent=>{
          const props=ent.properties;
          const stationType=props.type?.getValue()||'雷達';
          const status=props.status?.getValue()||'active';
          const hexCol=SSN_TYPE_COLORS[stationType]||'#00E5FF';
          const col=Cesium.Color.fromCssColorString(hexCol);
          const isRetired=(status==='decommissioned');
          ent.billboard=undefined;
          ent.point=new Cesium.PointGraphics({
            pixelSize:isRetired?5:9,
            color:isRetired?col.withAlpha(0.45):col,
            outlineColor:Cesium.Color.BLACK.withAlpha(0.7),
            outlineWidth:1,
          });
          // 設定標籤
          const nameVal=props.name?.getValue()||'';
          ent.label=new Cesium.LabelGraphics({
            text:nameVal,
            font:'10px Tahoma,sans-serif',
            fillColor:isRetired?Cesium.Color.fromCssColorString('#9E9E9E'):Cesium.Color.WHITE,
            outlineColor:Cesium.Color.BLACK,
            outlineWidth:2,
            style:Cesium.LabelStyle.FILL_AND_OUTLINE,
            pixelOffset:new Cesium.Cartesian2(0,-13),
            distanceDisplayCondition:new Cesium.DistanceDisplayCondition(0,isRetired?2e6:7e6),
            show:true,
          });
          // 設定 infoBox 描述
          const locVal=props.location?.getValue()||'';
          const notesVal=props.notes?.getValue()||'';
          const statusLabel=isRetired?'<span style="color:#f44336">已除役</span>':'<span style="color:#4caf50">運作中</span>';
          ent.description=new Cesium.ConstantProperty(
            '<div style="font-family:Tahoma,sans-serif;font-size:13px;padding:8px;background:#1e1e1e;color:#e8e8e8">'
            +'<p style="color:#8ab4f8;font-weight:bold;font-size:14px;margin:0 0 8px">'+nameVal+'</p>'
            +'<table style="border-collapse:collapse;width:100%;font-size:12px">'
            +'<tr><td style="color:#aaa;padding:2px 10px 2px 0;white-space:nowrap">類型</td>'
            +'<td style="color:#e8e8e8"><b>'+stationType+'</b></td></tr>'
            +'<tr><td style="color:#aaa;padding:2px 10px 2px 0">位置</td>'
            +'<td style="color:#e8e8e8">'+locVal+'</td></tr>'
            +'<tr><td style="color:#aaa;padding:2px 10px 2px 0">狀態</td>'
            +'<td>'+statusLabel+'</td></tr>'
            +(notesVal?'<tr><td style="color:#aaa;padding:2px 10px 2px 0;vertical-align:top">備註</td>'
            +'<td style="color:#ccc">'+notesVal+'</td></tr>':'')
            +'</table></div>'
          );
          // ── 指示線：抬高標記點 + 垂直細線，防止 zoom-in 時被地形吞沒 ──────
          const posCart=ent.position.getValue(Cesium.JulianDate.now());
          if(posCart){
            const carto=Cesium.Cartographic.fromCartesian(posCart);
            const lon=Cesium.Math.toDegrees(carto.longitude);
            const lat=Cesium.Math.toDegrees(carto.latitude);
            const h=20000; // 20 km above ellipsoid，足以超越任何地形高點
            // 重新定位標記與標籤至 20 km 高度
            ent.position=new Cesium.ConstantPositionProperty(
              Cesium.Cartesian3.fromDegrees(lon,lat,h));
            // 地面(0m) → 標記點(20km) 垂直指示線，加入同一 DataSource 方便整批移除
            ds.entities.add({
              polyline:{
                positions:Cesium.Cartesian3.fromDegreesArrayHeights([lon,lat,0,lon,lat,h]),
                width:1,
                material:new Cesium.ColorMaterialProperty(
                  (isRetired
                    ? Cesium.Color.fromCssColorString('#9E9E9E')
                    : col
                  ).withAlpha(isRetired?0.30:0.55)
                ),
                clampToGround:false,
                arcType:Cesium.ArcType.NONE,
                distanceDisplayCondition:new Cesium.DistanceDisplayCondition(0,7e6),
              }
            });
          }
        });
      }else if(type==='vatsim'){
        // VATSIM FIR/UIR 邊界
        ds=await Cesium.GeoJsonDataSource.load(ref.url,{describe:function(){ return ''; },clampToGround:false});
        [...ds.entities.values].forEach(ent=>{
          const props=ent.properties;
          const firId=(props&&props.id)?props.id.getValue()||'':'';
          const rawOce=(props&&props.oceanic!=null)?props.oceanic.getValue():0;
          const oceanic=(rawOce===1||rawOce===true||rawOce==='1');
          const hexCol=oceanic?'#4FC3F7':'#80DEEA';
          const col=Cesium.Color.fromCssColorString(hexCol).withAlpha(oceanic?0.55:0.80);
          if(ent.polygon){
            ent.polygon.fill=new Cesium.ConstantProperty(false);
            ent.polygon.outline=new Cesium.ConstantProperty(true);
            ent.polygon.outlineColor=new Cesium.ConstantProperty(col);
            ent.polygon.outlineWidth=new Cesium.ConstantProperty(oceanic?1.0:1.5);
            ent.polygon.arcType=new Cesium.ConstantProperty(Cesium.ArcType.GEODESIC);
          }
          if(ent.polyline){
            ent.polyline.material=new Cesium.ColorMaterialProperty(col);
            ent.polyline.width=new Cesium.ConstantProperty(oceanic?1.0:1.5);
            ent.polyline.arcType=new Cesium.ConstantProperty(Cesium.ArcType.GEODESIC);
          }
          // 使用 label_lon/label_lat 將標籤置於 FIR 中心
          const lblLon=(props&&props.label_lon!=null)?props.label_lon.getValue():null;
          const lblLat=(props&&props.label_lat!=null)?props.label_lat.getValue():null;
          if(firId){
            if(lblLon!=null&&lblLat!=null){
              ent.position=new Cesium.ConstantPositionProperty(
                Cesium.Cartesian3.fromDegrees(Number(lblLon),Number(lblLat),0));
            }
            ent.label=new Cesium.LabelGraphics({
              text:firId,
              font:'11px "Segoe UI",sans-serif',
              fillColor:Cesium.Color.fromCssColorString(hexCol),
              outlineColor:Cesium.Color.BLACK,
              outlineWidth:2,
              style:Cesium.LabelStyle.FILL_AND_OUTLINE,
              distanceDisplayCondition:new Cesium.DistanceDisplayCondition(0,3e6),
              pixelOffset:new Cesium.Cartesian2(0,0),
              show:true,
            });
          }
          ent.description=new Cesium.ConstantProperty(
            '<div style="font-family:Tahoma,sans-serif;font-size:13px;padding:8px;background:#1e1e1e;color:#e8e8e8">'
            +'<p style="color:'+hexCol+';font-weight:bold;font-size:14px;margin:0 0 8px">FIR: '+firId+'</p>'
            +'<table style="border-collapse:collapse;width:100%;font-size:12px">'
            +'<tr><td style="color:#aaa;padding:2px 10px 2px 0">類型</td>'
            +'<td style="color:#e8e8e8">'+(oceanic?'洋區 FIR (Oceanic)':'陸區 FIR (Continental)')+'</td></tr>'
            +(lblLat!=null&&lblLon!=null
              ?'<tr><td style="color:#aaa;padding:2px 10px 2px 0">標籤座標</td>'
              +'<td style="color:#e8e8e8">'+Number(lblLat).toFixed(2)+'&deg;, '+Number(lblLon).toFixed(2)+'&deg;</td></tr>':'')
            +'</table></div>'
          );
        });
      }
      await viewer.dataSources.add(ds);
      ref.set(ds);
    }catch(e){
      cb.checked=false;
      console.warn('向量圖層載入失敗: '+type, e);
    }
  }else{
    const ds=ref.get();
    if(ds) viewer.dataSources.remove(ds,true);
    ref.set(null);
  }
}

// ── 工具函數 ───────────────────────────────────────────────────────────────
function doRefresh(){ fetchPos(); }

function setUpdateInterval(ms){
  clearInterval(arTimer);
  arTimer=setInterval(fetchPos, ms);
}

// ── 啟動 ───────────────────────────────────────────────────────────────────
init().then(()=>{
  // 預設 60 秒自動更新（與 radio 預設值一致）
  arTimer=setInterval(fetchPos, 60000);
}).catch(e=>{
  document.getElementById('stat').textContent='初始化失敗: '+e.message;
  console.error(e);
});
</script>
</body>
</html>"""


# ─────────────────────────────────────────────────────────────────────────────
# 底圖目錄
# ─────────────────────────────────────────────────────────────────────────────
_BASEMAP_CATALOG: dict[str, dict] = {
    "default": {
        "file": "globe_texture.jpg", "label": "預設底圖",
        "credit": "NASA Blue Marble", "type": "single",
        "candidates": [_SCRIPT_DIR / "data" / "globe_texture.jpg"],
    },
    "2k": {
        "file": "land_shallow_topo_2048.jpg", "label": "藍色大理石 2K",
        "credit": "NASA Visible Earth", "type": "single",
        "candidates": [_TEXTURE_DIR / "land_shallow_topo_2048.jpg"],
    },
    "topo": {
        "file": "world_topo_bathy_5400.jpg", "label": "地形+海底地形 5.4K",
        "credit": "NASA Visible Earth", "type": "single",
        "candidates": [_TEXTURE_DIR / "world_topo_bathy_5400.jpg"],
    },
    "natural_earth": {
        "file": "", "label": "Natural Earth II",
        "credit": "Natural Earth / CesiumJS (已內建)", "type": "tms",
        "tms_url": "/cesium/Assets/Textures/NaturalEarthII/",
        "candidates": [
            _CESIUM_LOCAL_DIR / "Assets" / "Textures" / "NaturalEarthII" / "tilemapresource.xml"
        ],
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Flask 應用程式
# ─────────────────────────────────────────────────────────────────────────────
def create_app() -> Flask:
    app = Flask(__name__)
    CORS(app, resources={r"/api/*": {"origins": "*"}})

    @app.get("/cesium/<path:filename>")
    def cesium_static(filename: str):
        safe = (_CESIUM_LOCAL_DIR / filename).resolve()
        if not str(safe).startswith(str(_CESIUM_LOCAL_DIR.resolve())):
            return make_response("Forbidden", 403)
        if not safe.is_file():
            return make_response(f"Cesium asset not found: {filename}", 404)
        return send_from_directory(str(_CESIUM_LOCAL_DIR), filename)

    @app.get("/js/<path:filename>")
    def js_static(filename: str):
        safe = (_JS_DIR / filename).resolve()
        if not str(safe).startswith(str(_JS_DIR.resolve())):
            return make_response("Forbidden", 403)
        if not safe.is_file():
            return make_response(f"JS file not found: {filename}", 404)
        return send_from_directory(str(_JS_DIR), filename)

    @app.get("/api/globe_texture/<string:name>")
    def api_globe_texture(name: str):
        entry = _BASEMAP_CATALOG.get(name)
        if not entry:
            return make_response("Unknown basemap", 404)
        for p in entry.get("candidates", []):
            if p.is_file():
                return send_from_directory(str(p.parent), p.name)
        return make_response("Texture not found", 404)

    @app.get("/api/textures")
    def api_textures():
        out = []
        for key, entry in _BASEMAP_CATALOG.items():
            etype = entry.get("type", "single")
            available = any(p.is_file() for p in entry.get("candidates", []))
            item: dict[str, Any] = {
                "key": key,
                "label": entry["label"],
                "credit": entry.get("credit", ""),
                "type": etype,
                "available": available,
            }
            if etype == "tms":
                item["tms_url"] = entry.get("tms_url", "")
            out.append(item)
        return jsonify(out)

    @app.get("/")
    def index():
        return _HTML, 200, {"Content-Type": "text/html; charset=utf-8"}

    @app.get("/api/config")
    def api_config():
        return jsonify({"cesium_token": CESIUM_ION_TOKEN, "has_token": bool(CESIUM_ION_TOKEN)})

    # ── Himawari API ──────────────────────────────────────────────────────────

    @app.get("/api/himawari/latest")
    def api_himawari_latest():
        """
        向 NICT 查詢最新 Himawari 時間戳記，組出 Band 16 圖像 URL 並回傳。
        回傳: {timestamp, image_url, band, available}
        若 NICT 無回應則 available=false，前端應顯示 fallback。
        """
        ts = _fetch_himawari_latest()
        if not ts:
            return jsonify({
                "available": False,
                "band": "B16",
                "timestamp": None,
                "image_url": None,
                "error": "無法取得 Himawari 最新時間戳記",
            })
        image_url = _build_b16_url(ts)
        return jsonify({
            "available": True,
            "band": "B16",
            "timestamp": ts["timestamp"],
            "image_url": image_url,
        })

    @app.get("/api/himawari/frame")
    def api_himawari_frame():
        """
        抓取最新 Himawari Band 16 全圓盤 PNG，重投影至 WGS84 等經緯度後回傳。

        重投影方法（smythdesign.com/blog/georeferencing-himawari/）：
          GEOS +proj=geos +h=35785863 +lon_0=140.7 +sweep=y
          全圓盤邊界：±5,500,000 m（11,000 px × 1,000 m/px ÷ 2）
          輸出：1440×720 等經緯度 WGS84 PNG，圓盤外透明（alpha=0）

        結果快取 600 秒（對應 Himawari 10 分鐘更新週期）。
        """
        ts = _fetch_himawari_latest()
        if ts:
            ts_key = ts.get("timestamp", "")

            # ── 快取命中 ──
            if ts_key and ts_key in _hw_reproj_cache:
                logger.debug("Himawari B16 快取命中: %s", ts_key)
                png_data = _hw_reproj_cache[ts_key]
                flask_resp = make_response(png_data)
                flask_resp.headers["Content-Type"]          = "image/png"
                flask_resp.headers["Cache-Control"]         = "public, max-age=600"
                flask_resp.headers["X-Himawari-Timestamp"]  = ts_key
                flask_resp.headers["X-Himawari-Reprojected"] = "cached"
                flask_resp.headers["Access-Control-Allow-Origin"] = "*"
                return flask_resp

            # ── 下載原始 NICT PNG ──
            img_url = _build_b16_url(ts)
            logger.info("Himawari B16 下載: %s", img_url)
            try:
                resp = requests.get(
                    img_url,
                    timeout=HIMAWARI_REQUESTS_TIMEOUT,
                    headers={"User-Agent": "ATRDC-TLE-Tracker/1.0"},
                )
                resp.raise_for_status()
                raw_png = resp.content

                # ── 重投影至 WGS84 ──
                if _REPROJECT_OK:
                    logger.info("Himawari B16 重投影中（GEOS → WGS84）…")
                    png_data = _reproject_himawari_to_wgs84(raw_png)
                    reproj_label = "geos->wgs84"
                    logger.info("重投影完成，輸出 %d bytes", len(png_data))
                else:
                    png_data = raw_png
                    reproj_label = "raw(pyproj_missing)"
                    logger.warning("pyproj/PIL/scipy 未安裝，回傳原始 GEOS 影像")

                # 更新記憶體快取（保留最新一筆即可）
                _hw_reproj_cache.clear()
                if ts_key:
                    _hw_reproj_cache[ts_key] = png_data

                # 寫入本地備份檔，供離線 fallback 使用
                try:
                    _HIMAWARI_CACHE_DIR.mkdir(parents=True, exist_ok=True)
                    _HIMAWARI_LOCAL_FILE.write_bytes(png_data)
                    logger.info("Himawari 本地備份已更新: %s", _HIMAWARI_LOCAL_FILE)
                except Exception as _save_exc:
                    logger.warning("Himawari 本地備份寫入失敗: %s", _save_exc)

                flask_resp = make_response(png_data)
                flask_resp.headers["Content-Type"]           = "image/png"
                flask_resp.headers["Cache-Control"]          = "public, max-age=600"
                flask_resp.headers["X-Himawari-Timestamp"]   = ts_key
                flask_resp.headers["X-Himawari-Band"]        = "B16"
                flask_resp.headers["X-Himawari-Reprojected"] = reproj_label
                flask_resp.headers["Access-Control-Allow-Origin"] = "*"
                return flask_resp
            except Exception as exc:
                logger.warning("Himawari B16 圖像抓取/重投影失敗: %s", exc)
        else:
            logger.warning("Himawari 最新時間戳記無法取得")

        # Fallback 1：使用本地端備份檔
        if _HIMAWARI_LOCAL_FILE.is_file():
            try:
                local_png = _HIMAWARI_LOCAL_FILE.read_bytes()
                mtime = _HIMAWARI_LOCAL_FILE.stat().st_mtime
                local_ts = datetime.fromtimestamp(mtime, tz=timezone.utc).strftime("%Y%m%d_%H%M%S")
                logger.info("Himawari 使用本地備份檔 (mtime=%s)", local_ts)
                flask_resp = make_response(local_png)
                flask_resp.headers["Content-Type"]          = "image/png"
                flask_resp.headers["Cache-Control"]         = "public, max-age=300"
                flask_resp.headers["X-Himawari-Source"]     = "local-file"
                flask_resp.headers["X-Himawari-LocalTime"]  = local_ts
                flask_resp.headers["Access-Control-Allow-Origin"] = "*"
                return flask_resp
            except Exception as exc:
                logger.warning("Himawari 本地備份讀取失敗: %s", exc)

        # Fallback 2：1×1 透明 PNG（完全離線且無備份時）
        _TRANSPARENT_PNG = (
            b"\x89PNG\r\n\x1a\n"
            b"\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
            b"\x08\x06\x00\x00\x00\x1f\x15\xc4\x89"
            b"\x00\x00\x00\nIDATx\x9cc\x00\x01\x00\x00\x05\x00\x01\r\n-\xb4"
            b"\x00\x00\x00\x00IEND\xaeB`\x82"
        )
        flask_resp = make_response(_TRANSPARENT_PNG)
        flask_resp.headers["Content-Type"] = "image/png"
        flask_resp.headers["Cache-Control"] = "no-cache"
        flask_resp.headers["X-Himawari-Status"] = "unavailable"
        flask_resp.headers["Access-Control-Allow-Origin"] = "*"
        return flask_resp

    # ── PRC 衛星 API（與 scenario03.py 相同）─────────────────────────────────

    @app.get("/api/prc/categories")
    def api_categories():
        counts: dict[str, int] = {}
        for info in prc_catalog.values():
            cat = info["category"]
            counts[cat] = counts.get(cat, 0) + 1
        return jsonify([
            {"category": cat, "count": cnt, "color": CATEGORY_COLORS.get(cat, "#78909C")}
            for cat, cnt in sorted(counts.items(), key=lambda x: -x[1])
        ])

    @app.get("/api/prc/positions")
    def api_positions():
        cats_param = request.args.get("cats", "").strip()
        selected   = (
            {c.strip() for c in cats_param.split(",") if c.strip()}
            if cats_param else None
        )
        target_ids = [
            nid for nid, info in prc_catalog.items()
            if selected is None or info["category"] in selected
        ]
        if not target_ids:
            return jsonify({"count": 0, "satellites": [],
                            "timestamp": datetime.now(timezone.utc).isoformat()})
        tle_map = get_tle_cache()
        results = []
        for nid in target_ids:
            info = prc_catalog[nid]
            tle  = tle_map.get(nid)
            if tle is None:
                continue
            pos = propagate_to_now(tle["line1"], tle["line2"])
            if pos is None:
                continue
            results.append({
                "norad_id": nid,
                "name":     info["name"],
                "category": info["category"],
                "color":    CATEGORY_COLORS.get(info["category"], "#78909C"),
                "lat":      round(pos["lat"], 4),
                "lon":      round(pos["lon"], 4),
                "alt_km":   pos["alt_km"],
            })
        return jsonify({"count": len(results), "satellites": results,
                        "timestamp": datetime.now(timezone.utc).isoformat()})

    @app.get("/api/prc/orbit")
    def api_orbit():
        try:
            norad_id = int(request.args.get("norad_id", 0))
        except ValueError:
            return jsonify({"error": "norad_id 必須為整數"}), 400
        hours = float(request.args.get("hours", "2"))
        pts   = min(int(request.args.get("pts", "120")), 720)
        if norad_id not in prc_catalog:
            return jsonify({"error": f"NORAD {norad_id} 不在 PRC 目錄中"}), 404
        tle = get_tle_cache().get(norad_id)
        if tle is None:
            return jsonify({"error": f"NORAD {norad_id} 無 TLE 資料"}), 404
        positions = propagate_arc(tle["line1"], tle["line2"], hours=hours, pts=pts)
        return jsonify({"norad_id": norad_id, "name": prc_catalog[norad_id]["name"],
                        "hours": hours, "pts": len(positions), "positions": positions})

    @app.get("/api/prc/profile/<int:norad_id>")
    def api_profile(norad_id: int):
        """讀取 sat_profiles/<norad_id>_*.md，轉成含 CSS 的 HTML 供 Cesium infoBox 顯示。"""
        info = prc_catalog.get(norad_id)
        if not info:
            return jsonify({"error": f"NORAD {norad_id} 不在 PRC 目錄中"}), 404

        profile_path: Path = info["profile_path"]
        try:
            raw_md = profile_path.read_text(encoding="utf-8")
        except OSError as exc:
            logger.error("無法讀取 profile 檔案 %s: %s", profile_path, exc)
            return jsonify({"error": "profile 檔案讀取失敗"}), 500

        if _MD_OK:
            html_body = _md.markdown(raw_md, extensions=["tables", "sane_lists"])
        else:
            import html as _h
            html_body = f'<pre style="white-space:pre-wrap;font-size:12px">{_h.escape(raw_md)}</pre>'

        return jsonify({
            "norad_id": norad_id,
            "name":     info["name"],
            "html":     _INFOBOX_CSS + html_body,
        })

    # ── 向量圖層 API ─────────────────────────────────────────────────────────

    @app.get("/api/layers/borders")
    def api_layer_borders():
        """
        代理抓取 Natural Earth 110m 國界 GeoJSON 並快取。

        ★ 只保留 geometry + ADMIN（國名），其餘 94 個欄位全部丟棄。
          Natural Earth 原始 GeoJSON 每個 feature 有 ~94 個 properties，
          Cesium GeoJsonDataSource 嘗試 enumerate 時會拋出
          RangeError: Too many properties to enumerate，
          因此必須在後端先 strip 再回傳。
        """
        import json as _json
        global _borders_cache
        if _borders_cache is None:
            _NE_URL = (
                "https://raw.githubusercontent.com/nvkelso/natural-earth-vector"
                "/master/geojson/ne_110m_admin_0_countries.geojson"
            )
            try:
                r = requests.get(_NE_URL, timeout=30,
                                 headers={"User-Agent": "ATRDC-TLE-Tracker/1.0"})
                r.raise_for_status()
                raw = _json.loads(r.content)
                # ── 只保留國名，丟棄所有其他欄位 ──────────────────────────────
                for feat in raw.get("features", []):
                    feat["properties"] = {}   # 拋棄全部屬性，只保留 geometry
                _borders_cache = _json.dumps(raw, separators=(",", ":")).encode("utf-8")
                logger.info("國界 GeoJSON 已精簡快取，%d bytes（原始 %d bytes）",
                            len(_borders_cache), len(r.content))
            except Exception as exc:
                logger.warning("無法載入 Natural Earth 國界: %s", exc)
                empty = b'{"type":"FeatureCollection","features":[]}'
                resp = make_response(empty)
                resp.headers["Content-Type"] = "application/json; charset=utf-8"
                return resp
        resp = make_response(_borders_cache)
        resp.headers["Content-Type"] = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = "public, max-age=86400"
        return resp

    @app.get("/api/layers/ssn_stations")
    def api_layer_ssn_stations():
        """回傳 SSN 地面觀測站 GeoJSON（hardcoded，依查證更新版 MD 整理）。"""
        import json as _json
        resp = make_response(_json.dumps(_SSN_STATIONS_GEOJSON, ensure_ascii=False).encode("utf-8"))
        resp.headers["Content-Type"] = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = "public, max-age=3600"
        return resp

    @app.get("/api/layers/vatsim_boundaries")
    def api_layer_vatsim_boundaries():
        """
        代理 VATSpy FIR/UIR Boundaries GeoJSON 並快取。
        來源: https://github.com/vatsimnetwork/vatspy-data-project

        精簡 properties，只保留 id / oceanic / label_lon / label_lat，
        避免 Cesium GeoJsonDataSource enumerate 大量欄位的潛在問題。
        """
        import json as _json
        global _vatsim_cache
        if _vatsim_cache is None:
            _VATSIM_URL = (
                "https://raw.githubusercontent.com/vatsimnetwork"
                "/vatspy-data-project/master/Boundaries.geojson"
            )
            try:
                r = requests.get(_VATSIM_URL, timeout=30,
                                 headers={"User-Agent": "ATRDC-TLE-Tracker/1.0"})
                r.raise_for_status()
                raw = _json.loads(r.content)
                _KEEP = {"id", "oceanic", "label_lon", "label_lat"}
                for feat in raw.get("features", []):
                    p = feat.get("properties") or {}
                    feat["properties"] = {k: p[k] for k in _KEEP if k in p}
                _vatsim_cache = _json.dumps(raw, separators=(",", ":")).encode("utf-8")
                logger.info("VATSpy FIR Boundaries 已快取：%d bytes（原始 %d bytes）",
                            len(_vatsim_cache), len(r.content))
            except Exception as exc:
                logger.warning("無法載入 VATSpy Boundaries: %s", exc)
                empty = b'{"type":"FeatureCollection","features":[]}'
                resp = make_response(empty)
                resp.headers["Content-Type"] = "application/json; charset=utf-8"
                return resp
        resp = make_response(_vatsim_cache)
        resp.headers["Content-Type"] = "application/json; charset=utf-8"
        resp.headers["Cache-Control"] = "public, max-age=3600"
        return resp

    @app.get("/api/logo")
    def api_logo():
        logo_path = _SCRIPT_DIR / "Logo_ATRDC.png"
        if not logo_path.exists():
            return "", 404
        with open(logo_path, "rb") as f:
            data = f.read()
        resp = make_response(data)
        resp.headers["Content-Type"] = "image/png"
        resp.headers["Cache-Control"] = "public, max-age=86400"
        return resp

    @app.errorhandler(Exception)
    def handle_error(err):
        from werkzeug.exceptions import HTTPException
        if isinstance(err, HTTPException):
            return err
        logger.exception("發生未預期的例外錯誤")
        return jsonify({"error": str(err)}), 500

    return app


app = create_app()

if __name__ == "__main__":
    logger.info(
        "Scenario 05-6 (3D Cesium + 台灣視角 + Himawari B16 + 向量圖層 + Radio 更新間隔 + LOGO) 啟動中 — http://%s:%d",
        HOST, PORT,
    )
    app.run(host=HOST, port=PORT, debug=True)
