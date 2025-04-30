import torch
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import xlsxwriter

# 自动检测设备
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

def generate_patterns(num_neurons, num_patterns):
    """生成固定数量的二值模式 (±1)"""
    # 修改1：显式指定为float32类型
    return (torch.randint(0, 2, (num_patterns, num_neurons), 
                         device=device, dtype=torch.float32) * 2 - 1)

def compute_weights(patterns):
    """基于Hebbian规则计算权重矩阵"""
    # 修改2：确保使用浮点运算
    W = torch.matmul(patterns.T, patterns)
    W.fill_diagonal_(0)
    return W / patterns.size(1)

def fixed_flip(state, num_flips):
    """固定数量神经元翻转"""
    flip_indices = torch.randperm(state.size(0), device=device)[:num_flips]
    state[flip_indices] *= -1
    return state

def async_update(W, state, max_iter=100):
    """异步更新直到收敛"""
    for _ in range(max_iter):
        update_order = torch.randperm(state.size(0), device=device)
        new_state = state.clone()
        for idx in update_order:
            new_state[idx] = torch.sign(W[idx] @ new_state)
        if torch.equal(new_state, state):
            break
        state = new_state
    return state

def simulate_network(N, params):
    """仿真流程"""
    patterns = generate_patterns(N, params['num_patterns'])
    W = compute_weights(patterns)
    
    # 基准准确率
    h = torch.matmul(W, patterns.T).T
    base_acc = torch.mean((torch.sign(h) == patterns).all(dim=1).float()).item()
    
    # 扰动实验
    perturb_acc = []
    for _ in range(params['num_perturb']):
        perturbed = fixed_flip(patterns.clone(), params['num_flips'])
        h_perturbed = torch.matmul(W, perturbed.T).T
        perturb_acc.append(
            torch.mean((torch.sign(h_perturbed) == patterns).all(dim=1).float()).item()
        )
    
    return base_acc, np.var(perturb_acc)

def main():
    params = {
        'num_patterns': 40,#固定存储的模式（记忆）数量
        'num_flips': 5,#每次扰动时固定翻转的神经元数量
        'num_perturb': 100,#每个网络尺寸的扰动实验次数
        'num_trials': 10,#每个网络尺寸的独立重复实验次数
        'N_range': range(30, 1001, 20)#要测试的神经元数量范围
    }
    
    results = {'N': [], 'Accuracy': [], 'Variance': []}
    
    for N in tqdm(params['N_range']):
        if N < params['num_flips']:
            continue
            
        acc_trials = []# 初始化一个空列表，用于存储每次实验的基准准确率（accuracy）
        var_trials = []# 初始化一个空列表，用于存储每次实验的扰动方差（variance）
        
        for _ in range(params['num_trials']):
            base_acc, perturb_var = simulate_network(N, params)
            acc_trials.append(base_acc)
            var_trials.append(perturb_var)
        
        results['N'].append(N)
        results['Accuracy'].append(np.mean(acc_trials))
        results['Variance'].append(np.mean(var_trials))

    # 创建Excel文件
    workbook = xlsxwriter.Workbook('C:\\Users\\zty\\Desktop\\预测相变\\程序\\hopfield\\output_xlsxwriter.xlsx')
    worksheet = workbook.add_worksheet()

    # 写入数据
    worksheet.write_column('A1', results['N'])
    worksheet.write_column('B1', results['Accuracy'])
    worksheet.write_column('C1', results['Variance'])

    workbook.close()
    print("数据已写入output_xlsxwriter.xlsx")
    
    # 可视化
    fig, ax1 = plt.subplots(figsize=(12, 6))
    ax1.plot(results['N'], results['Accuracy'], 'b-o')
    ax1.set_xlabel('Number of Neurons')
    ax1.set_ylabel('Accuracy', color='b')
    ax2 = ax1.twinx()
    ax2.plot(results['N'], results['Variance'], 'r-s')
    ax2.set_ylabel('Variance', color='r')
    plt.title("Hopfield Network Stability Analysis")
    plt.show()

if __name__ == "__main__":
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    main()