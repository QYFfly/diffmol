import os

def split_sdf(input_sdf_path, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    count = 0
    sdf_files = [f for f in os.listdir(input_sdf_path) if f.endswith('.sdf')]
    sdf_files.sort()  # 确保文件按顺序处理

    for sdf_file in sdf_files:
        with open(os.path.join(input_sdf_path, sdf_file), 'r') as file:
            sdf_content = file.read()
        
        sdf_blocks = sdf_content.split('$$$$\n')
        
        for idx, sdf_block in enumerate(sdf_blocks):
            if sdf_block.strip():  # 确保块不为空
                output_file_path = f"{output_dir}/{count+1}.sdf"
                with open(output_file_path, 'w') as output_file:
                    output_file.write(sdf_block + '$$$$\n')
                count += 1
    
    print(f"Total SDF files processed: {count}")

input_sdf_path = '/data/qinyifei/model/test_mol/smol/1/'  # 替换为你的输入SDF文件路径
output_dir = '/data/qinyifei/model/test_mol/smol/GEOM-Drugs/sdf/'  # 替换为你的输出目录
# input_sdf_path = '/data/qinyifei/model/out/JODO/test/'  # 替换为你的输入SDF文件路径
# output_dir = '/data/qinyifei/model/out/JODO/test/'  # 替换为你的输出目录
split_sdf(input_sdf_path, output_dir)
