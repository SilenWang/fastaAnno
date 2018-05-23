# Readme

## 脚本简介
因为客户需要写的将二代测序获得的变异展示在参考基因组序列上的小脚本, 目前功能:
- 将注释好的变异(snp / indel)标注在参考基因组上
- 将没有覆盖(depth 0)的位点标注出来

## 优化方向
- 将变异输入文件变更为vcf/avinput
- 变更未覆盖的输入文件方式?
- 对indel采取悬浮展示?
- 如果输入vcf是否可读入vcf的详细信息然后也以悬浮窗进行展示

## 目前结果展示
![demo_pic](https://github.com/SilenWang/fastaAnno/blob/master/pic/demo.jpeg)