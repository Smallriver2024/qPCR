这是一个基于 **Streamlit** 的网页工具，用于快速计算 qPCR 实验的相对表达量（ΔΔCt 方法）。
作者：
bilibili@ 外科小小硕
---

## 功能特色
- **数据输入灵活**：支持复制粘贴表格（Tab/逗号分隔）或上传 CSV/Excel 文件
- <img width="396" height="438" alt="image" src="https://github.com/user-attachments/assets/7544bc93-0fba-4f8a-8c6c-7390e63ac142" />
数据整理成这样子，分组那一列的名称一定要是group，每一行是一个样本，然后后面可以放内参和各种基因，同时也支持多基因导入。
整理好后直接复制到左边的格子里面，也可以直接上传文件。
- **对照组可配置**：可手动指定对照组（如 NC / Control）  
- **自动识别基因列**：可选择目标基因 Ct 列与管家基因 Ct 列  
- **自动计算**：  
  - ΔCt = 目标基因Ct − 管家基因Ct  
  - ΔΔCt = ΔCt − 对照组均值(ΔCt)  
  - 相对表达量 = 2^(-ΔΔCt)  
- **结果展示**：  
  - 每个样本的计算结果  
  - 分组均值 ± SD（可选）  
- **结果导出**：一键下载 Excel 文件，包含结果明细、参数设置和分组统计
- <img width="1320" height="752" alt="4241756114109_ pic" src="https://github.com/user-attachments/assets/08245280-4321-4716-9ec2-874bb8ff4eae" />
