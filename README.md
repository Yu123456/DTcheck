# DTcheck
DTcheck 是一个基于爬虫查询 KEGG, DrugBank, ChEMBL  3 个数据库验证药靶作用关系的小工具，需要输入药物靶标对的 KEGG DRUG ID, KEGG GENE ID。
针对 Yamanishi Y, et al 于 2008 年发布的 drug-target 数据集中 unknown-interaction 对，返回药靶作用对是否已经被上述 3 个数据库验证（即是否
已经可以从该数据库中查询到作用关系）。

### 论文信息
余冬华,郭茂祖,刘晓燕. 药物靶标作用关系预测结果评价及查询验证研究.

### 代码及数据
在 DTcheck/Data/ 文件夹下有需要输入的原始数据，代码在 DTcheck/Code/ 文件夹下

**注意** drugTargetCheck.py 中下述 4 条命令，需要针对待查询数据集而注释掉其余 3 个数据集，直接运行 drugTargetCheck.py 
 dn = dataName[1]      # GPCR 数据集
 dn = dataName[2]      # IC 数据集
 dn = dataName[0]      # Enzyme 数据集
 dn = dataName[3]      # NR 数据集
 
 ### E-mail: donghuayu@hit.edu.cn or guomaozu@bucea.edu.cn
