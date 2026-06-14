---
title: RhynoW SSA Dashboard
emoji: 🛰️
colorFrom: blue
colorTo: indigo
sdk: docker
pinned: false
app_port: 7860
---

# 太空態勢儀表板 / Space Situational Awareness Dashboard

**離線 SSA/SDA 視覺化儀表板**，以 Flask + CesiumJS 1.114 建構，無需外部網路即可運行。

![Python 3.11+](https://img.shields.io/badge/Python-3.11%2B-blue)
![Flask 3](https://img.shields.io/badge/Flask-3.x-lightgrey)
![CesiumJS 1.114](https://img.shields.io/badge/CesiumJS-1.114-green)
![Port 7860](https://img.shields.io/badge/Port-7860-orange)

---

## 功能概覽

### 主儀表板（3D 地球，`/`）

| 功能 | 說明 |
|------|------|
| **全目錄分類顯示** | 依國家 / 用途 / 年代 / 星座篩選，點選分類後顯示全量衛星位置（無抽樣） |
| **向量化 SGP4** | 使用 `sgp4.SatrecArray` 批次傳播，Starlink 全量（~6000 顆）< 1 s |
| **KD-tree 近距離掃描** | ECI 空間建樹，找出當前瞬間距離 ≤ 閾值的配對 |
| **物件搜尋** | NORAD ID 精確比對 / 名稱子串搜尋，top-20 即時顯示並定位 |
| **全球國界** | Natural Earth 110m GeoJSON，離線可用 |
| **SSN 觀測站** | 17 個 Space Surveillance Network 地面站，含類型 / 狀態 / 備註 |
| **完全離線** | CesiumJS 從本地 `data/cesium/` 載入；地球貼圖本地供圖 |

### 台北覆蓋分析（2D 地圖，`/taipei`）

| 功能 | 說明 |
|------|------|
| **台北 2000 km 覆蓋圈** | Cesium 2D 穩定顯示，無飄移問題 |
| **四類衛星分類顯示** | 美國商用光學、中國商用光學、中國軍用偵察、台灣 TASA |
| **仰角 / 方位計算** | ECI → ECEF → ENU 向量化計算，標示可見衛星（仰角 > 5°） |
| **24 小時過頂預報** | SGP4 逐分鐘批次傳播，找出升弧 / 降弧 / 最大仰角時刻 |
| **±30 天時間軸** | 滑桿回溯過去（讀取 DB 最近 TLE）/ 預測未來（SGP4 外推） |
| **歷史 / 現在 / 預測 模式** | 自動標示當前時間模式，「現在」模式每 60 秒自動更新 |
| **衛星名稱標籤** | Zoom ≥ 7 時自動顯示衛星名稱（`DistanceDisplayCondition`） |

---

## 快速開始

### 步驟 1 — 環境安裝

```bash
git clone https://github.com/<your-org>/scenario04-advanced01
cd scenario04-advanced01

python -m venv .venv
.venv\Scripts\activate        # Windows
# source .venv/bin/activate   # Linux / macOS

pip install -r requirements.txt
```

### 步驟 2 — 下載 CesiumJS（一次性）

前往 [CesiumJS Releases](https://github.com/CesiumGS/cesium/releases/tag/1.114) 下載
`Cesium-1.114.zip`，解壓後將 `Build/Cesium/` 目錄內容放置於：

```
data/cesium/          ← 此目錄內應含 Cesium.js、Widgets/、Workers/ 等
```

驗證路徑是否正確：

```
data/cesium/Cesium.js           ✓
data/cesium/Widgets/widgets.css ✓
data/cesium/Workers/            ✓
```

### 步驟 3 — 準備資料庫

將 TLE DuckDB 資料庫（`space_db_slim.duckdb` 或 `space_db.duckdb`）放置於專案根目錄：

```
scenario04-advanced01/
└── space_db_slim.duckdb    ← 放在這裡
```

資料庫需包含以下資料表：

| 資料表 | 說明 |
|--------|------|
| `raw_tle_archive` | TLE 歷史記錄（含 `norad_id`, `object_name`, `line1`, `line2`, `epoch_utc`） |
| `sat_n2yo_metadata` | 衛星元資料（`source_code`, `launch_date`, `intl_code`） |

> 如需從 Space-Track 下載並建立資料庫，請參考
> [ATRDC_TLE_Tracker](https://github.com/<your-org>/ATRDC_TLE_Tracker)。

### 步驟 4 — 設定環境變數

```bash
copy .env.example .env    # Windows
# cp .env.example .env    # Linux / macOS
```

| 變數 | 預設值 | 說明 |
|------|--------|------|
| `DB_PATH` | `./space_db_slim.duckdb` | DuckDB 資料庫路徑 |
| `CESIUM_ION_TOKEN` | *(空)* | Cesium ION token（[免費申請](https://ion.cesium.com/tokens)，不設定仍可運行） |
| `HOST` | `0.0.0.0` | Flask 監聽位址 |
| `PORT` | `5013` | Flask 監聽埠 |
| `CONJ_THRESHOLD_KM` | `10.0` | 近距離掃描閾值（km） |
| `CONJ_TTL` | `120` | 接近事件快取存活秒數 |
| `STATS_TTL` | `600` | 統計快取存活秒數 |
| `LOG_LEVEL` | `INFO` | 日誌等級 |

### 步驟 5 — 啟動

```bash
python scenario04-Cesium-advanced02.py
```

開啟瀏覽器：

- 主儀表板（3D）：**http://localhost:7860**
- 台北覆蓋分析（2D）：**http://localhost:7860/taipei**

---

## 目錄結構

```
scenario04-advanced01/
├── scenario04-Cesium-advanced02.py   # Flask 應用（後端 + 前端 HTML 單檔）
├── sat_metadata.csv                  # 衛星元資料補充（名稱 / 國家 / 星座覆蓋）
├── requirements.txt
├── .env.example
└── data/
    ├── borders.geojson               # Natural Earth 110m 國界（已隨附）
    ├── globe_texture.jpg             # NASA Blue Marble 地球貼圖（已隨附）
    └── cesium/                       # CesiumJS 1.114（git-ignored，手動下載）
        └── .gitkeep
```

---

## API 端點

### 主儀表板

| 方法 | 路徑 | 說明 |
|------|------|------|
| GET | `/api/stats` | 統計（總數、國家、用途、年代、星座分布） |
| GET | `/api/positions` | `?ftype=country&fval=美國` → 批次傳播並回傳位置 |
| GET | `/api/position/<norad_id>` | 單顆衛星當前位置 |
| GET | `/api/conjunctions` | `?threshold_km=10&max_pairs=200` → KD-tree 近距離掃描 |
| GET | `/api/search` | `?q=<NORAD ID 或名稱>` → top-20 搜尋結果 |
| GET | `/api/globe_texture` | 地球貼圖（JPEG）|
| GET | `/api/layers/borders` | 國界 GeoJSON |
| GET | `/api/layers/ssn_stations` | SSN 觀測站 GeoJSON |
| GET | `/cesium/<path>` | 本機 CesiumJS 靜態資源 |

### 台北覆蓋分析

| 方法 | 路徑 | 說明 |
|------|------|------|
| GET | `/taipei` | Cesium 2D 台北覆蓋頁面 |
| GET | `/api/taipei_coverage` | 當前四類衛星位置 + 仰角 / 方位 |
| GET | `/api/taipei_coverage_at` | `?ts=<ISO 時間>` → 指定時刻覆蓋狀況 |
| GET | `/api/taipei_passes` | 未來 24 小時過頂預報（四類） |
| GET | `/api/taipei_passes_at` | `?ts=<ISO 時間>` → 指定時刻起算的過頂預報 |

---

## 衛星分類（台北覆蓋頁面）

| 類別 ID | 名稱 | 關鍵字 |
|---------|------|--------|
| `US_EO` | 美國商用光學衛星 | WorldView、GeoEye、Legion、SkySat、Pelican |
| `CN_COMM` | 中國商用光學衛星 | SuperView、吉林、珠海、高分 |
| `CN_MIL` | 中國軍用偵察衛星 | 遙感 Yaogan、尖兵 Jianbing |
| `TW_TASA` | 台灣 TASA 衛星 | Formosat-5/-7/-8、COSMIC-2 |

---

## 架構說明

### Cesium 2D 穩定化

台北覆蓋頁使用以下設定避免畫面向下飄移：

- 直接以 `SCENE2D` 初始化（不走 morph 轉場）
- `camera.setView()`（同步）而非 `flyTo()`（非同步動畫）
- `viewer.clock.shouldAnimate = false`
- 停用 `screenSpaceCameraController` 的 `enableRotate` / `enableTilt`

### 並發控制

前端使用以下機制防止頻繁操作造成 Cesium 渲染錯誤：

- `AbortController` 取消過期的 fetch 請求
- `_loading` flag 防止並發 `_loadForTs()` 呼叫
- 按鈕 1 秒冷卻時間（`_lastAction`）
- 滑桿 500 ms debounce
- `satDs.entities.suspendEvents()` / `resumeEvents()` 批次實體更新

### 向量化 SGP4

使用 `sgp4.SatrecArray` 批次傳播（需 sgp4 ≥ 2.0）；傳入前以
`np.ascontiguousarray()` 確保記憶體連續性：

```python
jds = np.ascontiguousarray(jd_fr[:, 0])
frs = np.ascontiguousarray(jd_fr[:, 1])
e, r, _ = SatrecArray([...]).sgp4(jds, frs)
```

---

## 相依套件

- **[sgp4](https://pypi.org/project/sgp4/)** ≥ 2.22 — 衛星傳播（含 `SatrecArray` 批次模式）
- **[duckdb](https://duckdb.org/)** ≥ 0.10 — 嵌入式分析資料庫（TLE 歷史查詢）
- **[scipy](https://scipy.org/)** ≥ 1.11 — `cKDTree` 空間索引
- **[numpy](https://numpy.org/)** ≥ 1.26 — 向量化座標轉換
- **[flask](https://flask.palletsprojects.com/)** ≥ 3.0 — REST API 框架
- **[CesiumJS](https://cesium.com/platform/cesiumjs/)** 1.114 — 3D/2D 地球視覺化

---

## 注意事項

- 近距離配對為「當前瞬間 ECI 位置距離」，並非完整的 TCA/Pc 計算。
  詳細接近分析（Stage A/B/C）請參考 `conjunction_pipeline.py`。
- SGP4 傳播精度隨 TLE 年齡降低；建議使用 7 天內的 TLE。
- 時間軸回溯功能需要 DB 內有歷史 TLE（`raw_tle_archive.epoch_utc`）；
  未來預測則以當前 TLE 進行 SGP4 外推，精度隨時間增加而下降。
- `SatrecArray` 批次模式需要 `sgp4 >= 2.0`；版本較舊時自動退回逐顆傳播。
