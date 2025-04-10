import json
from pathlib import Path
from typing import Dict, Any
from utils import (
    validity_connect,
    calculate_uniqueness,
    computely_novelty,
    process_unique_sdf,
    calculate_QED_for_subfolders,
    calculate_SA_for_subfolders,
    calculate_overall_diversity,
    calculate_scaffold_diversity,
    calculate_structural_properties_2D,
    extract_structural_properties,
    calculate_stable,
    calculate_FCD,
    calculate_OPLS3_preparation,
    calculate_MMFF_preparation,
    calculate_bond_angles,
    calculate_dihedral_angles,
    calculate_bond_length,
    calculate_JSD_bond_length,
    calculate_JSD_angles,
    calculate_JSD_dihedral_angles,
    MolecularMetrics
)

# 常量定义
# 获取当前脚本所在目录
SCRIPT_DIR = Path(__file__).parent.absolute()

# 使用相对路径定位config.json
CONFIG_PATH = SCRIPT_DIR / "config.json"
DEFAULT_FORCEFIELD = "MMFF"
DATASET_NAME = "Drug"
PROJECT_ROOT = Path(__file__).parent.absolute()
class DrugMetricsCalculator:
    """Drug指标计算器"""
    
    def __init__(self, input_dir: Path, output_dir: Path):
        """初始化"""
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.config = self._load_config()
        self.forcefield = self.config.get("forcefield", DEFAULT_FORCEFIELD)
        self.ligprep_executable = self.config.get("ligprep_executable")
        # 设置默认路径
        # self.default_train_smiles = Path("/data/qinyifei/model/Diff-3D/data/data_reference/Drug/Drug_smiles.txt")
        # self.default_training_set = Path("/data/qinyifei/model/Diff-3D/data/training_set_reference/Drug/smiles.smi")
        # self.default_fpscores =Path("/data/qinyifei/model/Diff-3D/data/fpscores.pkl.gz")
        self.default_train_smiles = PROJECT_ROOT / "data/data_reference/Drug/Drug_smiles.txt"
        self.default_training_set = PROJECT_ROOT / "data/training_set_reference/Drug/smiles.smi"
        self.default_fpscores = PROJECT_ROOT / "data/fpscores.pkl.gz"
        # 从config获取或使用默认路径
        config_fpscores= self.config.get("fpscores", "")
        self.fpscores = config_fpscores if config_fpscores else self.default_fpscores

        config_smiles_path = self.config.get("train_smiles_file", "")
        self.train_smiles_file = Path(config_smiles_path) if config_smiles_path else self.default_train_smiles
        
        config_training_set = self.config.get("training_set_file", "")
        self.training_set_file = Path(config_training_set) if config_training_set else self.default_training_set

        self.calc_qed = self.config.get("calculate_qed", True)
        self.calc_sa = self.config.get("calculate_sa", True)

        self.calc_2d_props = self.config.get("calculate_2d_properties", True)
        self.calc_stability = self.config.get("calculate_stability", True)
        self.calc_fcd = self.config.get("calculate_fcd", True)
        
        self.calc_overall_div = self.config.get("calculate_overall_diversity", True)
        self.calc_scaffold_div = self.config.get("calculate_scaffold_diversity", True)
        
        self.calc_geometry = self.config.get("calculate_geometry_metrics", True)
        
        # 初始化分子指标计算器
        self.molecular_metrics = MolecularMetrics()
    
    def _load_config(self) -> Dict[str, Any]:
        """加载配置文件"""
        try:
            with open(CONFIG_PATH, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            raise FileNotFoundError(f"Config file not found at {CONFIG_PATH}")
        except json.JSONDecodeError:
            raise ValueError(f"Invalid JSON format in config file {CONFIG_PATH}")
    
    def _validate_directories(self) -> None:
        """验证输入输出目录"""
        if not self.input_dir.exists():
            raise FileNotFoundError(f"Input directory not found: {self.input_dir}")
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def _check_opls3_requirements(self) -> None:
        """检查OPLS3所需配置"""
        if self.forcefield == "OPLS3" and not self.config.get("ligprep_executable"):
            raise ValueError("OPLS3 requires ligprep_executable path in config")
    
    def calculate_basic_metrics(self) -> None:
        """计算基础指标"""
        print("Calculating basic metrics...")
        print(f"Using training smiles file: {self.train_smiles_file}")
        
        # 分子有效性计算
        self.molecular_metrics.process_and_save(self.input_dir, self.output_dir)
        
        # 其他基础指标
        calculate_uniqueness(self.output_dir, DATASET_NAME)
        computely_novelty(self.output_dir, DATASET_NAME, str(self.train_smiles_file))
        process_unique_sdf(self.input_dir, self.output_dir, DATASET_NAME)
    
    def calculate_chemical_properties(self) -> None:
        """计算化学属性（根据配置决定是否执行）"""
        print("\nCalculating chemical properties...")
        
        # QED计算
        if self.calc_qed:
            print("Calculating QED...")
            calculate_QED_for_subfolders(self.output_dir, DATASET_NAME)
        else:
            print("Skipping QED calculation (disabled in config)")
        
        # SA计算
        if self.calc_sa:
            print("Calculating SA...")
            calculate_SA_for_subfolders(self.output_dir, DATASET_NAME, self.fpscores)
        else:
            print("Skipping SA calculation (disabled in config)")
    
    
    def calculate_diversity_metrics(self) -> None:
            """计算多样性指标"""
            print("\nCalculating diversity metrics...")
            if self.calc_overall_div:
                print("Calculating overall diversity...")
                calculate_overall_diversity(self.output_dir, DATASET_NAME)
            else:
                print("Skipping overall diversity (disabled in config)")
            
            if self.calc_scaffold_div:
                print("Calculating scaffold diversity...")
                calculate_scaffold_diversity( self.output_dir, DATASET_NAME)
            else:
                print("Skipping scaffold diversity (disabled in config)")
    
    def calculate_2d_properties(self) -> None:
        """计算2D属性(根据配置决定是否执行)"""
        print("\nCalculating 2D properties...")
        
        # 2D结构计算（共享一个开关）
        if self.calc_2d_props:
            print("Calculating structural properties...")
            calculate_structural_properties_2D(self.output_dir, DATASET_NAME)
            extract_structural_properties(self.output_dir, DATASET_NAME)
        else:
            print("Skipping 2D structural properties (disabled in config)")
        
        # 稳定性计算
        if self.calc_stability:
            print("Calculating stability...")
            calculate_stable(self.output_dir, DATASET_NAME)
        else:
            print("Skipping stability calculation (disabled in config)")
        
        # FCD计算
        if self.calc_fcd:
            print("Calculating FCD...")
            calculate_FCD(str(self.training_set_file), self.output_dir, DATASET_NAME)
        else:
            print("Skipping FCD calculation (disabled in config)")
        
    def calculate_geometry_metrics(self) -> None:
            """计算所有几何相关指标（统一开关控制）"""
            if not self.calc_geometry:
                print("Skipping all geometry-related metrics (disabled in config)")
                return
            
            print("\n===== Calculating Geometry Metrics =====")
            
            # 力场相关计算
            print(f"\n[1/3] Forcefield ({self.forcefield}) preparation...")
            if self.forcefield == "OPLS3":
                calculate_OPLS3_preparation(self.output_dir, DATASET_NAME, self.ligprep_executable)
            else:
                calculate_MMFF_preparation(self.output_dir, DATASET_NAME)
            
            # 几何参数计算
            print("\n[2/3] Geometric parameters...")
            calculate_bond_length(self.output_dir, DATASET_NAME, self.forcefield)
            calculate_bond_angles(self.output_dir, DATASET_NAME, self.forcefield)
            calculate_dihedral_angles(self.output_dir, DATASET_NAME, self.forcefield)
            
            # 相似性指标
            print("\n[3/3] Similarity metrics...")
            calculate_JSD_bond_length(self.output_dir, DATASET_NAME)
            calculate_JSD_angles(self.output_dir, DATASET_NAME)
            calculate_JSD_dihedral_angles(self.output_dir, DATASET_NAME)
            
            print("\nAll geometry metrics completed.")
    
    def run_all_metrics(self) -> None:
        """运行所有指标计算"""
        try:
            self._validate_directories()
            self._check_opls3_requirements()
            
            print("\n" + "="*50)
            print(f"Starting {DATASET_NAME} metrics calculation with config:")
            print(f"- Forcefield: {self.forcefield}")
            print(f"- Training SMILES: {self.train_smiles_file}")
            print(f"- Training Set: {self.training_set_file}")
            if self.forcefield == "OPLS3":
                print(f"- Ligprep executable: {self.config.get('ligprep_executable', 'NOT SET')}")
            print("="*50 + "\n")
            
            self.calculate_basic_metrics()
            self.calculate_chemical_properties()
            self.calculate_diversity_metrics()
            self.calculate_2d_properties()
            self.calculate_geometry_metrics()  # 替换原来的两个方法
            
            print(f"{DATASET_NAME} metrics processing completed successfully.")
        except Exception as e:
            print(f"Error in {DATASET_NAME} metrics processing: {str(e)}")
            raise

def drug_metrics(input_directory: str, output_directory: str) -> None:
    """Drug指标计算入口函数"""
    calculator = DrugMetricsCalculator(
        input_dir=Path(input_directory),
        output_dir=Path(output_directory))
    calculator.run_all_metrics()