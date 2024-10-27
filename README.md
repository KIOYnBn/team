#  readme

## 库的声明

* python == 3.9.19
* pytorch = 2.4.1
* pandas == 2.2.2
* numpy == 1.26.4
* bio == 1.6.2
* icecream为个人调试使用，不影响运行
## 声明
* 数据存贮格式以DUD-Z数据集为标准数据集
* data_acquire.py旨在从pdb获取序列信息与三维结构信息
  * 该库独立于模型之外，属于数据收集部分，仅运行一次即可。

## 工作日志

* 26, 10, 24
  * 完成data_acquire的BIO序列获取
