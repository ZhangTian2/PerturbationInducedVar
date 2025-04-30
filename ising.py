import torch
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import xlsxwriter
import os

class IsingModel2D:
    def __init__(self, size=50, h0=0.0, J=1.0, device='cuda' if torch.cuda.is_available() else 'cpu'):
        """
        初始化2D Ising模型

        参数:
            size: 晶格尺寸 (size x size)
            h0: 恒定外磁场
            J: 最近邻相互作用强度
            device: 'cuda' 或 'cpu'
        """
        self.size = size
        self.h0 = h0
        self.J = J
        self.device = device

        # 初始化自旋矩阵
        self.spins = 2 * torch.randint(0, 2, (size, size), device=device) - 1

        # 最近邻卷积核
        self.neighbor_kernel = torch.tensor([[0, 1, 0],
                                             [1, 0, 1],
                                             [0, 1, 0]], dtype=torch.float32, device=device)

    def energy(self, h_noise=0.0):
        """
        计算系统总能量
        参数:
            h_noise: 噪声磁场（默认 0）
        """
        neighbor_sum = torch.nn.functional.conv2d(
            self.spins.unsqueeze(0).unsqueeze(0).float(),
            self.neighbor_kernel.unsqueeze(0).unsqueeze(0),
            padding=1
        ).squeeze()

        interaction_energy = -0.5 * self.J * torch.sum(self.spins * neighbor_sum)
        field_energy = -(self.h0 + h_noise) * torch.sum(self.spins)

        return interaction_energy + field_energy

    def magnetization(self):
        return torch.mean(self.spins.float())

    def metropolis_step(self, temperature, h_noise=0.0):
        num_flips = self.size * self.size
        flip_indices = torch.randint(0, self.size, (2, num_flips), device=self.device)

        neighbor_sum = torch.nn.functional.conv2d(
            self.spins.unsqueeze(0).unsqueeze(0).float(),
            self.neighbor_kernel.unsqueeze(0).unsqueeze(0),
            padding=1
        ).squeeze()

        selected_spins = self.spins[flip_indices[0], flip_indices[1]]
        delta_E = 2 * selected_spins * (
            self.J * neighbor_sum[flip_indices[0], flip_indices[1]] +
            self.h0 + h_noise
        )

        accept_prob = torch.exp(-delta_E / temperature)
        accept = torch.rand(num_flips, device=self.device) < accept_prob

        self.spins[flip_indices[0, accept], flip_indices[1, accept]] *= -1

    def simulate(self, temperature, steps=1000, equil_steps=500, h_noise=0.0):
        for _ in range(equil_steps):
            self.metropolis_step(temperature, h_noise)

        mag_sum = 0.0
        for _ in range(steps - equil_steps):
            self.metropolis_step(temperature, h_noise)
            mag_sum += abs(self.magnetization())

        return mag_sum / (steps - equil_steps)

def simulate_with_noise_analysis(
    size=50, h0=0.0, J=1.0, ha_std=0.01,
    min_temp=0.5, max_temp=4.0, temp_steps=30,
    mc_steps=1000, samples_per_temp=10):

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")

    temperatures = torch.linspace(min_temp, max_temp, temp_steps, device=device)
    magnetizations_no_noise = torch.zeros(temp_steps, device=device)
    magnetization_var_with_noise = torch.zeros(temp_steps, device=device)

    for i, temp in enumerate(tqdm(temperatures, desc="Simulating")):
        # 模型1：无噪声，计算平均磁化强度
        model_clean = IsingModel2D(size=size, h0=h0, J=J, device=device)
        magnetizations_no_noise[i] = model_clean.simulate(temp, steps=mc_steps)

        # 模型2：多次采样有噪声，计算方差
        mags_with_noise = []
        for _ in range(samples_per_temp):
            h_noise = float(np.random.normal(loc=0.0, scale=ha_std))
            model_noisy = IsingModel2D(size=size, h0=h0, J=J, device=device)
            mag = model_noisy.simulate(temp, steps=mc_steps, h_noise=h_noise)
            mags_with_noise.append(mag.item())

        magnetization_var_with_noise[i] = torch.tensor(mags_with_noise, device=device).var(unbiased=True)

    # 转移到CPU
    temps_np = temperatures.cpu().numpy()
    mags_np = magnetizations_no_noise.cpu().numpy()
    vars_np = magnetization_var_with_noise.cpu().numpy()

    # 创建Excel文件
    workbook = xlsxwriter.Workbook('C:\\Users\\zty\\Desktop\\预测相变\\程序\\Ising模型\\output_xlsxwriter.xlsx')
    worksheet = workbook.add_worksheet()

    # 写入数据
    worksheet.write_column('A1', temps_np)
    worksheet.write_column('B1', mags_np)
    worksheet.write_column('C1', vars_np)

    workbook.close()
    print("数据已写入output_xlsxwriter.xlsx")

    # 绘图
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    plt.plot(temps_np, mags_np, 'o-', label='Avg Magnetization (no noise)')
    plt.xlabel("Temperature")
    plt.ylabel("Magnetization")
    plt.title("Average Magnetization vs Temperature")
    plt.grid(True)
    plt.axvline(x=2.269, color='r', linestyle='--', label='Critical Temp ~2.27')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(temps_np, vars_np, 'o-', color='orange', label='Magnetization Variance (with noise)')
    plt.xlabel("Temperature")
    plt.ylabel("Variance")
    plt.title("Magnetization Variance due to Noise")
    plt.grid(True)
    plt.axvline(x=2.269, color='r', linestyle='--')
    plt.legend()

    plt.tight_layout()
    plot_path = os.path.join('C:\\Users\\zty\\Desktop\\预测相变\\程序\\Ising模型', 'ising_perturbation_analysis.png')
    plt.savefig(plot_path)
    plt.show()
    print(f"图表已保存到 {plot_path}")

# 运行主函数
if __name__ == "__main__":
    simulate_with_noise_analysis(
        size=50,
        h0=1.0,           # 恒定磁场
        J=5.0,            # 相互作用强度
        ha_std=0.01,      # 噪声磁场的标准差
        min_temp=0.5,
        max_temp=40,
        temp_steps=50,
        mc_steps=1000,
        samples_per_temp=10
    )
