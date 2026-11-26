"""Open a local native directory dialog outside Streamlit's script thread."""

from pathlib import Path
import subprocess
import sys


def pick_folder(current):
    result = subprocess.run(
        [sys.executable, '-X', 'utf8', str(Path(__file__).resolve()), current],
        capture_output=True, encoding='utf-8', check=True,
        creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0),
    )
    return result.stdout.strip()


def main(current):
    from tkinter import Tk, filedialog

    start = Path(current).expanduser()
    if not start.is_dir():
        start = Path.home()
    window = Tk()
    try:
        window.withdraw()
        window.attributes('-topmost', True)
        selected = filedialog.askdirectory(parent=window, initialdir=str(start),
                                            mustexist=True, title='选择本地仿真结果目录')
    finally:
        window.destroy()
    print(str(Path(selected)) if selected else '')


if __name__ == '__main__':
    main(sys.argv[1])
