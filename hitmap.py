import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.animation import FuncAnimation, PillowWriter
import os


def make_stream_video(stream_csv):
    """시간 순서(EventID)에 따라 hit가 누적되는 GIF 생성"""
    if not os.path.exists(stream_csv):
        return

    try:
        df = pd.read_csv(stream_csv)
        if df.empty: return

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlim(15.5, -0.5)
        ax.set_ylim(-0.5, 15.5)

        ax.set_xlabel("Column ID (15 -> 0)")
        ax.set_ylabel("Row ID (15 <- 0)")
        ax.grid(True, linestyle='--', alpha=0.4)

        scat = ax.scatter([], [], s=180, c=[], cmap='plasma', vmin=0, vmax=25, edgecolors='black', alpha=0.7)
        xs, ys, colors = [], [], []

        def update(frame):
            xs.append(df.iloc[frame]["Col"])
            ys.append(df.iloc[frame]["Row"])
            colors.append(df.iloc[frame]["TOA_ns"])

            scat.set_offsets(list(zip(xs, ys)))
            scat.set_array(colors)
            ax.set_title(
                f"Hit Stream Evolution | EventID = {int(df.iloc[frame]['EventID'])} | Current TOA = {colors[-1]:.2f} ns")
            return scat,

        ani = FuncAnimation(fig, update, frames=len(df), interval=150, blit=True)
        base, _ = os.path.splitext(stream_csv)
        out_name = base + "_stream.gif"

        ani.save(out_name, writer=PillowWriter(fps=6))
        print(f"-> Saved Stream GIF: {out_name}")
    except Exception as e:
        print(f"[Warning] Failed to encode stream GIF for {stream_csv}: {e}")
    finally:
        plt.close()


def plot_global_heatmap(heatmap_csv):
    """개별 채널(Port) 16x16 누적 Occupancy 히트맵 시각화"""
    if not os.path.exists(heatmap_csv): return

    try:
        df = pd.read_csv(heatmap_csv, header=None)
        if df.empty or df.shape[1] != 16: return

        # 매트릭스 방향 정렬 (오른쪽 아래가 0,0)
        df = df.iloc[::-1].reset_index(drop=True)
        df = df.iloc[:, ::-1]

        plt.figure(figsize=(9, 7))
        sns.heatmap(df, cmap='viridis', annot=False, cbar_kws={'label': 'Total Hits'},
                    yticklabels=list(range(15, -1, -1)), xticklabels=list(range(15, -1, -1)))

        plt.yticks(rotation=0)
        plt.title(f"Global Occupancy Map (Bottom-Right = 0,0)\n({os.path.basename(heatmap_csv)})")
        plt.xlabel("Column ID (15 -> 0)")
        plt.ylabel("Row ID (15 -> 0)")

        base, _ = os.path.splitext(heatmap_csv)
        out_name = base + "_plot.png"
        plt.savefig(out_name, dpi=300)
        plt.close()
        print(f"-> Saved Individual Heatmap: {out_name}")
    except Exception as e:
        pass


def plot_combined_heatmaps(dir_path):
    """
    Port_A와 Port_B의 히트맵 데이터를 더하여 하나의 16x16 그리드에 합산 플롯 생성
    (우측 하단이 Col 0, Row 0)
    """
    # [수정] 파일명 탐색을 CH_A/CH_B 에서 Port_A/Port_B 로 변경
    port_a_file = os.path.join(dir_path, "global_Port_A_heatmap.csv")
    port_b_file = os.path.join(dir_path, "global_Port_B_heatmap.csv")

    if not (os.path.exists(port_a_file) and os.path.exists(port_b_file)):
        return

    try:
        # 1. 두 파일 읽어오기
        df_a = pd.read_csv(port_a_file, header=None)
        df_b = pd.read_csv(port_b_file, header=None)

        if df_a.shape[1] != 16 or df_b.shape[1] != 16:
            return

        # 2. 두 채널(Port)의 히트수 행렬 더하기 (Pandas DataFrame 덧셈)
        df_combined = df_a + df_b

        # 3. 매트릭스 방향 정렬 (오른쪽 아래가 0,0)
        # 행(Row) 상하 반전 후, 열(Col) 좌우 반전
        df_combined = df_combined.iloc[::-1].reset_index(drop=True).iloc[:, ::-1]

        # 4. 하나의 그림으로 시각화
        plt.figure(figsize=(9, 7))

        sns.heatmap(
            df_combined,
            cmap='viridis',
            annot=False,
            cbar_kws={'label': 'Total Hits (Port_A + Port_B)'}, # [수정] 라벨 이름 변경
            yticklabels=list(range(15, -1, -1)),
            xticklabels=list(range(15, -1, -1))
        )

        plt.yticks(rotation=0)

        # [수정] 타이틀 명칭 변경
        plt.title(f"Combined Global Occupancy Map (Port_L + Port_R)\nBottom-Right = 0,0")
        plt.xlabel("Column ID (15 -> 0)")
        plt.ylabel("Row ID (15 -> 0)")

        out_name = os.path.join(dir_path, "global_Combined_Heatmap_plot.png")
        plt.savefig(out_name, dpi=300)
        plt.close()
        print(f"-> Saved Single Combined Heatmap: {out_name}")

    except Exception as e:
        print(f"[Error] Failed to generate combined heatmap: {e}")

if __name__ == "__main__":
    print("=== ETROC Python Visualization Tool ===")
    target_dir = "results"

    if not os.path.exists(target_dir):
        print(f"[Error] Target folder '{target_dir}' is missing. Run C++ parser first.")
        exit(1)

    # 1. stream gif (첫 번째 파일 하나만)
    stream_files = sorted([f for f in os.listdir(target_dir) if "hit_stream" in f and f.endswith(".csv")])
    if len(stream_files) > 0:
        sample_file = os.path.join(target_dir, stream_files[0])
        make_stream_video(sample_file)

    # 2. 개별 global heatmap
    heatmap_files = sorted(
        [f for f in os.listdir(target_dir) if "global" in f and "heatmap" in f and f.endswith(".csv")])
    for hm_file in heatmap_files:
        plot_global_heatmap(os.path.join(target_dir, hm_file))

    # 3. Port_A & Port_B 통합 global heatmap
    plot_combined_heatmaps(target_dir)

    print("=======================================")