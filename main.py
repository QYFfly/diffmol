import os
os.environ["OMP_NUM_THREADS"] = "1"
import json
from pathlib import Path
from qm9_metrics import qm9_metrics  # 从 qm9_metrics.py 导入主函数
from drug_metrics import drug_metrics  # 从 drug_metrics.py 导入主函数

# 加载配置文件
def load_config(config_path):
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config

# 主函数
def main():
    # 读取配置文件路径
    # 获取当前脚本所在目录
    SCRIPT_DIR = Path(__file__).parent.absolute()

    # 使用相对路径定位config.json
    config_path= SCRIPT_DIR / "config.json"
    # config_path = '/data/qinyifei/model/Diff-3D/config.json'  # 配置文件路径
    config = load_config(config_path)  # 加载配置

    # 获取 dataset、input_directory 和 output_directory
    dataset = config.get('dataset', 'QM9')  # 默认使用 'QM9' 数据集
    input_directory = config.get('input_directory', '/data/qinyifei/model/Diff-3D/data/')
    output_directory = config.get('output_directory', '/data/qinyifei/model/Diff-3D/output/')

    # 根据 dataset 来调用不同的处理模块
    if dataset == 'QM9':
        print("Using QM9 metrics...")
        qm9_metrics(input_directory, output_directory)  # 调用 qm9_metrics.py 的处理函数
    elif dataset == 'Drug':
        print("Using Drug metrics...")
        drug_metrics(input_directory, output_directory)  # 调用 drug_metrics.py 的处理函数
    else:
        print(f"Unknown dataset: {dataset}. Please specify either 'QM9' or 'Drug'.")

if __name__ == "__main__":
    main()
