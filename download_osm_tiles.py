#!/usr/bin/env python3
"""
download_osm_tiles.py
---------------------
下載 OpenStreetMap 圖磚 z0-8，儲存為 MBTiles (SQLite) 格式。

來源：   tile.openstreetmap.org
授權：   © OpenStreetMap contributors (ODbL)
政策：   https://operations.osmfoundation.org/policies/tiles/
         本腳本用於教育/離線個人使用，附 Attribution，限速符合 ToS。

輸出：   data/tiles/osm_z0-8.mbtiles  (~700 MB)
用法：   python download_osm_tiles.py [--zoom 0-8] [--workers 8] [--resume]
"""

from __future__ import annotations

import argparse
import sqlite3
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

# ── 設定 ──────────────────────────────────────────────────────────────────────
ZOOM_MIN   = 0
ZOOM_MAX   = 8
OUTPUT     = Path(__file__).parent / "data" / "tiles" / "osm_z0-8.mbtiles"
TILE_URL   = "https://tile.openstreetmap.org/{z}/{x}/{y}.png"
USER_AGENT = (
    "ATRDC-TLE-Tracker/1.0 offline-tile-download "
    "(educational use; +https://github.com/ATRDC)"
)
MAX_WORKERS   = 8     # 並發下載數（OSM ToS：合理限速）
REQUEST_DELAY = 0.05  # 每張下載後最短等待 (s)，8 worker × 20/s ≈ 160 req/s max
BATCH_COMMIT  = 200   # 每幾張 commit 一次
RETRY_MAX     = 3


# ── MBTiles 初始化 ─────────────────────────────────────────────────────────────
def init_db(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(str(path), check_same_thread=False)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("""CREATE TABLE IF NOT EXISTS metadata
                    (name TEXT PRIMARY KEY, value TEXT)""")
    conn.execute("""CREATE TABLE IF NOT EXISTS tiles (
                    zoom_level  INTEGER NOT NULL,
                    tile_column INTEGER NOT NULL,
                    tile_row    INTEGER NOT NULL,
                    tile_data   BLOB    NOT NULL,
                    PRIMARY KEY (zoom_level, tile_column, tile_row))""")
    meta = [
        ("name",        "OpenStreetMap z0-8"),
        ("format",      "png"),
        ("minzoom",     str(ZOOM_MIN)),
        ("maxzoom",     str(ZOOM_MAX)),
        ("description", "OpenStreetMap offline tiles for ATRDC TLE Tracker"),
        ("attribution", "© OpenStreetMap contributors"),
    ]
    conn.executemany("INSERT OR REPLACE INTO metadata VALUES (?,?)", meta)
    conn.commit()
    return conn


def existing_tiles(conn: sqlite3.Connection) -> set[tuple[int,int,int]]:
    """已下載的 (z, x, tms_y) 集合，用於 resume。"""
    rows = conn.execute(
        "SELECT zoom_level, tile_column, tile_row FROM tiles"
    ).fetchall()
    return set(rows)


# ── 下載邏輯 ──────────────────────────────────────────────────────────────────
def download_tile(
    session: requests.Session,
    z: int, x: int, y: int,
) -> bytes | None:
    url = TILE_URL.format(z=z, x=x, y=y)
    for attempt in range(RETRY_MAX):
        try:
            resp = session.get(url, timeout=20)
            if resp.status_code == 200:
                return resp.content
            if resp.status_code == 429:
                time.sleep(5 * (attempt + 1))
        except requests.RequestException:
            time.sleep(2 * (attempt + 1))
    return None


# ── 主程式 ────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(description="Download OSM tiles to MBTiles")
    parser.add_argument("--zoom",    default=f"{ZOOM_MIN}-{ZOOM_MAX}")
    parser.add_argument("--workers", type=int, default=MAX_WORKERS)
    parser.add_argument("--output",  default=str(OUTPUT))
    parser.add_argument("--resume",  action="store_true", default=True,
                        help="跳過已下載的圖磚（預設啟用）")
    args = parser.parse_args()

    z_min, z_max = map(int, args.zoom.split("-"))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # 計算總數
    total = sum(4**z for z in range(z_min, z_max + 1))
    print(f"目標：z{z_min}-z{z_max}，共 {total:,} 張圖磚")
    print(f"輸出：{out_path}")
    print(f"並發：{args.workers} workers")

    conn = init_db(out_path)
    lock = threading.Lock()

    # Resume：取得已完成清單
    done_set: set[tuple[int,int,int]] = set()
    if args.resume:
        done_set = existing_tiles(conn)
        if done_set:
            print(f"Resume 模式：已有 {len(done_set):,} 張，跳過。")

    # 建立任務清單（XYZ 坐標，MBTiles 存 TMS y = 2^z-1 - y）
    tasks: list[tuple[int,int,int]] = []
    for z in range(z_min, z_max + 1):
        n = 2**z
        for x in range(n):
            for y in range(n):
                tms_y = n - 1 - y
                if (z, x, tms_y) not in done_set:
                    tasks.append((z, x, y))

    remaining = len(tasks)
    already   = total - remaining
    print(f"待下載：{remaining:,} 張（已完成 {already:,} 張）")
    if remaining == 0:
        print("全部已完成！")
        conn.close()
        return

    # 下載
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    done_count = already
    fail_count = 0
    batch: list[tuple] = []
    t0 = time.monotonic()

    def worker(z: int, x: int, y: int) -> tuple[int,int,int,bytes|None]:
        data = download_tile(session, z, x, y)
        time.sleep(REQUEST_DELAY)
        return z, x, y, data

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(worker, z, x, y): (z,x,y) for z,x,y in tasks}
        for fut in as_completed(futures):
            z, x, y, data = fut.result()
            if data:
                tms_y = (2**z - 1) - y
                batch.append((z, x, tms_y, data))
                done_count += 1
            else:
                fail_count += 1

            if len(batch) >= BATCH_COMMIT:
                with lock:
                    conn.executemany(
                        "INSERT OR REPLACE INTO tiles VALUES (?,?,?,?)", batch
                    )
                    conn.commit()
                    batch.clear()

            # 進度
            processed = done_count + fail_count - already
            if processed % 500 == 0 or processed == remaining:
                elapsed = time.monotonic() - t0
                rate    = processed / elapsed if elapsed > 0 else 0
                eta_min = (remaining - processed) / rate / 60 if rate > 0 else 0
                pct     = 100 * done_count / total
                print(
                    f"  [{pct:5.1f}%] {done_count:,}/{total:,} 完成"
                    f"  失敗 {fail_count}"
                    f"  {rate:.1f} t/s"
                    f"  ETA {eta_min:.0f} min"
                )

    # 最後批次 commit
    if batch:
        conn.executemany(
            "INSERT OR REPLACE INTO tiles VALUES (?,?,?,?)", batch
        )
        conn.commit()

    conn.close()
    size_mb = out_path.stat().st_size / 1024**2
    elapsed = time.monotonic() - t0
    print(f"\n完成！{done_count:,}/{total:,} 張，失敗 {fail_count} 張")
    print(f"檔案大小：{size_mb:.1f} MB，耗時 {elapsed/60:.1f} 分鐘")
    print(f"輸出：{out_path}")


if __name__ == "__main__":
    main()
