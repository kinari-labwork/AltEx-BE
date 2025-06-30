import argparse

def main():
    parser = argparse.ArgumentParser(
        description="テストです",
    )

    # 明示的に -v/--version を追加
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='0.1.0',
        help='バージョン情報を表示して終了します'
    )

    parser.parse_args()

    print("実行中...（特に何もしません）")

if __name__ == '__main__':
    main()