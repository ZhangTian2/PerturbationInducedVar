import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader
from torchvision import datasets, transforms
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score
import xlsxwriter
import os

# 检查CUDA是否可用
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# 加载MNIST数据集
transform = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize((0.1307,), (0.3081,))
])

train_dataset = datasets.MNIST(root='./data', train=True, download=True, transform=transform)
test_dataset = datasets.MNIST(root='./data', train=False, download=True, transform=transform)

# 定义数据加载器
train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=1000, shuffle=False)

# 定义RBM类
class RBM(nn.Module):
    def __init__(self, visible_dim, hidden_dim):
        super(RBM, self).__init__()
        self.W = nn.Parameter(torch.randn(hidden_dim, visible_dim, device=device) * 0.1)
        self.h_bias = nn.Parameter(torch.zeros(hidden_dim, device=device))
        self.v_bias = nn.Parameter(torch.zeros(visible_dim, device=device))
    
    def sample_from_p(self, p):
        return torch.relu(torch.sign(p - torch.rand_like(p)))
    
    def v_to_h(self, v):
        p_h = torch.sigmoid(F.linear(v, self.W, self.h_bias))
        return p_h, self.sample_from_p(p_h)
    
    def h_to_v(self, h):
        p_v = torch.sigmoid(F.linear(h, self.W.t(), self.v_bias))
        return p_v, self.sample_from_p(p_v)
    
    def forward(self, v, k=1):
        h0, h_sample = self.v_to_h(v)
        for _ in range(k):
            v_recon, v_sample = self.h_to_v(h_sample)
            h_recon, h_sample = self.v_to_h(v_sample)
        return v, v_recon
    
    def free_energy(self, v):
        vbias_term = v.mv(self.v_bias)
        wx_b = F.linear(v, self.W, self.h_bias)#v @ self.W.t() + self.h_bias
        hidden_term = wx_b.exp().add(1).log().sum(1)
        return (-hidden_term - vbias_term).mean()#自由能的完整数学表达式：F v  = -∑ᵢ aᵢvᵢ - ∑ⱼ ln(1 + exp(∑ᵢ Wⱼᵢvᵢ + bⱼ))

# 定义分类模型
class RBMClassifier(nn.Module):
    def __init__(self, visible_dim, hidden_dim, num_classes):
        super(RBMClassifier, self).__init__()
        self.rbm = RBM(visible_dim, hidden_dim)
        self.classifier = nn.Linear(hidden_dim, num_classes, device=device)#该层实现的变换为：output=x@W+b, 其中：x：输入（形状为[batch_size, hidden_dim]的隐藏层表示）W：可学习权重矩阵（形状为[num_classes, hidden_dim]）b：可学习偏置（形状为[num_classes]）
    
    def perturb_hidden(self, h, num_flips, n_perturbations):
        batch_size, hidden_dim = h.size()
        flip_indices = torch.randint(0, hidden_dim, 
                                    (n_perturbations, batch_size, num_flips),
                                    device=device)
        flip_mask = torch.zeros(n_perturbations, batch_size, hidden_dim,
                              dtype=torch.bool, device=device)
        flip_mask.scatter_(2, flip_indices, True)
        h_expanded = h.unsqueeze(0).expand(n_perturbations, -1, -1)
        h_perturbed = h_expanded.clone()
        h_perturbed[flip_mask] = 1 - h_perturbed[flip_mask]
        return h_perturbed.reshape(-1, hidden_dim)
    
    def forward(self, x, perturb_flips=0, n_perturbations=1):
        x = x.view(-1, 28*28)
        _, h = self.rbm.v_to_h(x)
        if perturb_flips > 0 and n_perturbations > 1:
            h = self.perturb_hidden(h, perturb_flips, n_perturbations)
        return self.classifier(h)

# 评估函数
def evaluate(model):
    model.eval()
    y_true, y_pred = [], []
    with torch.no_grad():
        for data, target in test_loader:
            data = data.to(device)
            data = (data > 0.5).float()  # 二值化处理
            output = model(data.view(-1, 28*28))
            pred = output.argmax(dim=1)
            y_true.extend(target.cpu().numpy())
            y_pred.extend(pred.cpu().numpy())
    return accuracy_score(y_true, y_pred)

# 扰动评估函数
def evaluate_perturbation(model, num_flips=5, n_perturbations=20):
    model.eval()
    all_preds = []
    y_true = []
    
    with torch.no_grad():
        for data, target in test_loader:
            data, target = data.to(device), target.to(device)
            data = (data > 0.5).float()  # 二值化处理
            batch_size = data.shape[0]
            
            # 生成扰动后的预测
            logits = model(data, perturb_flips=num_flips, n_perturbations=n_perturbations)
            logits = logits.view(n_perturbations, batch_size, -1)
            preds = logits.argmax(dim=2).cpu().numpy()
            
            # 收集真实标签和预测结果
            all_preds.append(preds)
            y_true.extend(target.cpu().numpy())  # 原始标签，不是重复的
    
    # 计算每个扰动版本的准确率
    all_preds = np.concatenate(all_preds, axis=1)
    y_true = np.array(y_true)
    accuracies = []
    
    # 对每个扰动版本计算准确率
    for i in range(n_perturbations):
        acc = accuracy_score(y_true, all_preds[i])
        accuracies.append(acc)
    
    return np.mean(accuracies), np.var(accuracies)

# 训练和评估函数
def train_and_evaluate(hidden_units, num_flips=5, n_perturbations=20):
    model = RBMClassifier(28*28, hidden_units, 10).to(device)
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.CrossEntropyLoss()
    
    # 预训练RBM
    print(f"\nPre-training RBM with {hidden_units} hidden units...")
    for epoch in range(5):
        for data, _ in tqdm(train_loader, desc=f"Pre-train Epoch {epoch+1}"):
            data = data.view(-1, 28*28).to(device)
            data = (data > 0.5).float()
            v, v_recon = model.rbm(data)
            loss = model.rbm.free_energy(v) - model.rbm.free_energy(v_recon)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
    
    # 训练分类器
    print(f"\nTraining classifier with {hidden_units} hidden units...")
    for epoch in range(10):
        model.train()
        for data, target in tqdm(train_loader, desc=f"Train Epoch {epoch+1}"):
            data, target = data.to(device), target.to(device)
            data = (data > 0.5).float()
            optimizer.zero_grad()
            output = model(data)
            loss = criterion(output, target)
            loss.backward()
            optimizer.step()
    
    # 评估
    original_acc = evaluate(model)
    perturb_mean, perturb_var = evaluate_perturbation(model, num_flips, n_perturbations)
    print(f"Original Acc: {original_acc:.4f}, Perturbed Mean: {perturb_mean:.4f}, Variance: {perturb_var:.6f}")
    
    return original_acc, perturb_mean, perturb_var

# 实验参数
# hidden_units_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34]

hidden_units_list = list(range(1, 102, 5))
num_flips = 5
n_perturbations = 20
# num_flips = 1
# n_perturbations = 50

# 运行实验
results = {
    'hidden_units': [],
    'original_acc': [],
    'perturb_mean': [],
    'perturb_var': []
}

for units in hidden_units_list:
    orig_acc, pm, pv = train_and_evaluate(units, num_flips, n_perturbations)
    results['hidden_units'].append(units)
    results['original_acc'].append(orig_acc)
    results['perturb_mean'].append(pm)
    results['perturb_var'].append(pv)


# 创建Excel文件
workbook = xlsxwriter.Workbook('C:\\Users\\zty\\Desktop\\output_xlsxwriter.xlsx')
worksheet = workbook.add_worksheet()

# 写入数据
worksheet.write_column('A1', results['hidden_units'],)
worksheet.write_column('B1', results['original_acc'])
worksheet.write_column('C1', results['perturb_mean'])
worksheet.write_column('D1', results['perturb_var'])

workbook.close()
print("数据已写入output_xlsxwriter.xlsx")

# 可视化
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(results['hidden_units'], results['original_acc'], 'o-', label='Original')
plt.plot(results['hidden_units'], results['perturb_mean'], 's--', label='Perturbed Mean')
plt.xscale('log')
plt.xlabel('Hidden Units')
plt.ylabel('Accuracy')
plt.title('Accuracy vs Network Size')
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(results['hidden_units'], results['perturb_var'], 'd-.', color='red')
plt.xscale('log')
plt.xlabel('Hidden Units')
plt.ylabel('Variance')
plt.title('Variance of Accuracy under Perturbations')
plt.grid(True)

plt.tight_layout()
plot_path = os.path.join('C:\\Users\\zty\\Desktop', 'rbm_perturbation_analysis.png')
plt.savefig(plot_path)
plt.show()
print(f"图表已保存到 {plot_path}")