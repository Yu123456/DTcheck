# coding:utf-8
# 2018-06-22
# spider
# 通过 drug ID 核对其对应的 target ID
# drug 需要在 3 个数据库中查找核对
# KEGG DRUG, ChEMBL, DrugBank

# 通过 target ID 核对其对应的 drug ID
# target 在 KEGG GENES 数据库中查找核对


import pandas as pd
import numpy as np
import requests
import re
import time

from functions import *



# 数据集名称
dataName = ['Enzyme', 'GPCR', 'IC','NR']
# dn = dataName[1]      # GPCR 数据集
# dn = dataName[2]      # IC 数据集
# dn = dataName[0]      # Enzyme 数据集
dn = dataName[3]      # NR 数据集

print('%s 处理中……' %dn)

outName1 = dn+'KEGGCheck.csv'
outName2 = dn+'KEGGUncheck.csv'
outName3 = dn+'DrugBankCheck.csv'
outName4 = dn+'DrugBankUncheck.csv'
outName5 = dn+'ChEMBLCheck.csv'
outName6 = dn+'ChEMBLUncheck.csv'
outName = dn+'ResultCheck.csv'     # 整理后的结果(只针对有已验证结果)
outName7 = dn+'TopChecked.csv'     # 在原始数据中添加查询结果
outName8 = dn+'hsaCheck.csv'       # 通过 target ID 查找到的已验证作用关系
outName9 = dn+'hsaUncheck.csv'     # target ID 读取网页失败的 ID
# header = (Drug ID, Target ID, score, label, ChEMBL, DrugBank, KEGG, TRUE)
# ChEMBL : C, DrugBank : D, KEGG : K
# 只要在任何一个数据库得到验证， TRUE : T

data = pd.read_csv(dn+'Top.csv')
# 由于读取进去的 columns name 有空格，先去掉
data.columns = data.columns.str.strip()

data['Target ID'] = data['Target ID'].str.strip()    # 去除 Target ID 列中 hsa id 的头尾空格
data['Drug ID'] = data['Drug ID'].str.strip()     # 去除 Drug ID 列中 drug id 的头尾空格

drugID = data['Drug ID']     # 提取 Drug ID 列数据
#  drug ID 去除重复并递增排序
drugID = sorted(drugID)
drugID = np.unique(drugID)


unCheckKEGG = []        # KEGG 读取网页失败的 drug id
CheckedKEGG = []        # KEGG 已经验证的 (drug id, target id), tuple 类型
unCheckDB = []         # DrugBank 网页读取失败的 drug id
CheckedDB = []         # DrugBank 已经验证的 (drug id, target id), tuple 类型
unCheckChEMBL = []     # ChEMBL 网页读取失败的 drug id
CheckedChEMBL = []     # ChEMBL 已经验证的 (drug id, target id), tuple 类型

# 为了调试代码，先运行前两个 Drug ID
# drugID = drugID[0:3]
# print 'drugID'
# print drugID

for id in drugID:
    print('Drug %s' %id)
    # 提取 drug id 所对应的所有 target id
    # Series 格式
    targetID = data.loc[data['Drug ID'] == id,'Target ID']
    # print 'targetID type'
    # print type(targetID)
    id_number = id[1:]

    # KEGG DRUG 数据库
    temp_hsa_kegg = []     # 在 KEGG 中找到对应 drug id 的 hsa id, list 格式
    # 爬虫请求头信息
    user_agent = 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:49.0) Gecko/20100101 Firefox/49.0'
    referer = 'www.kegg.jp'
    headers = {'User-Agent': user_agent, 'Referer': referer}
    # 请求链接
    url = 'https://www.kegg.jp/dbget-bin/www_bget?dr:D{}'.format(id_number)
    # 获取网页
    flag, response = get_page(url=url, headers=headers)
    # 如果网页获取成功，就进行 HSA ID 提取，否则将其加入 unCheckKEGG
    if flag:
        # 获取网页成功
        print '处理 KEGG   数据库'
        # 1、Target 栏
        # 调用函数，解析 Target 栏
        target_id = get_target_id(response=response.text)
        if target_id:
            # target_id 非空
            temp_hsa_kegg.extend(target_id)

        # 2、Pathway 栏
        # 调用函数，解析 Pathway 栏
        pathway_id = get_pathway_id(response=response.text)
        if pathway_id:
            # pathway_id 非空
            temp_hsa_kegg.extend(pathway_id)

        # 3、 Metabolism 栏
        # 调用函数，解析 Metabolism 栏
        metabolism_id = get_metabolism_id(response=response.text)
        if metabolism_id:
            # metabolism_id 非空
            temp_hsa_kegg.extend(metabolism_id)

        # 4、 Interaction 栏
        # 调用函数，解析 Interaction 栏
        interaction_id = get_interaction_id(response=response.text)
        if interaction_id:
            # interaction_id 非空
            temp_hsa_kegg.extend(interaction_id)

        # 去除重复
        if temp_hsa_kegg:
            # temp_hsa_kegg 非空，可进行去重处理
            temp_hsa_kegg = np.unique(temp_hsa_kegg)  # 在 KEGG DRUG 中查询到的 drug 所对应的所有 hsa id

        # print 'temp_hsa_kegg'
        # print temp_hsa_kegg

        # 核查药靶作用关系
        drugtarget_check = check_interaction(drugid=id,hsa_predict=targetID,hsa_db=temp_hsa_kegg)
        if drugtarget_check:
            CheckedKEGG.extend(drugtarget_check)         # 将查找到的作用关系保存到 CheckedKEGG 中

        # print 'CheckedKEGG'
        # print CheckedKEGG

        # 5、获取 ChEMBL, DrugBank 中的两个链接
        chembl_url, drugbank_url = get_url(response.text)

        print('处理 DrugBank 数据库')
        # DrugBank 数据库
        temp_hsa_db = []  # 在 DrugBank 中找到对应 drug id 的 hsa id, list 格式
        # 爬虫请求头信息
        user_agent = 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:49.0) Gecko/20100101 Firefox/49.0'
        referer = 'www.drugbank.ca'
        headers = {'User-Agent': user_agent, 'Referer': referer}
        if drugbank_url:
            # drugbank url 非空（由于只有一个 drug id, 无需循环读取链接）
            url = drugbank_url[0]
            # 获取网页
            flag, response = get_page(url=url, headers=headers)
            if flag:
                # 网页获取成功
                uniprot_id = get_uniprot_id(response.text)
                if uniprot_id:
                    # uniprot_id 非空
                    # 说明存在对应的 uniprot id, 去除重复
                    uniprot_id = np.unique(uniprot_id)
                    # 对查询到的 uniprot id 转换成 KEGG id (hsa id)
                    temp_hsa_db = IDmapping(u_id=uniprot_id)
                    # print('temp_hsa_db')
                    # print temp_hsa_db
                    # 核查药靶作用关系
                    drugtarget_check = check_interaction(drugid=id, hsa_predict=targetID, hsa_db=temp_hsa_db)
                    if drugtarget_check:
                        CheckedDB.extend(drugtarget_check)  # 将查找到的作用关系保存到 CheckedDB 中

            else:
                # DrugBank 网页获取失败
                print('DrugBank 网页获取失败')
                unCheckDB.append(id)

        print('处理 ChEMBL 数据库')
        # ChEMBL 数据库
        temp_hsa_chembl = []  # 在 ChEMBL 中找到对应 drug id 的 hsa id, list 格式

        # 在所有的 chembl_url 中的链接，也可能指向相同的 hsa id
        for url in chembl_url:
            # 逐条链接处理
            temp_id = get_hsaid_chembl(url=url)
            if temp_id:
                # temp_id 非空
                temp_hsa_chembl.extend(temp_id)

        drugtarget_check = check_interaction(drugid=id, hsa_predict=targetID, hsa_db=temp_hsa_chembl)
        if drugtarget_check:
            CheckedChEMBL.extend(drugtarget_check)  # 将查找到的作用关系保存到 CheckedChEMBL 中

        # print('temp_hsa_chembl')
        # print temp_hsa_chembl

    else:
        # 获取网页失败
        print('KEGG 网页获取失败')
        unCheckKEGG.append(id)


# 通过 Target ID 验证 Drug ID
targetID = data['Target ID']  # 提取 Target ID 列数据
# target ID 去除重复并递增排序
targetID = sorted(targetID)
targetID = np.unique(targetID)

unCheckKEGGhsa = []    # KEGG 读取网页失败的 target id
CheckedKEGGhsa = []    # KEGG 已验证的 (drug id, target id), tuple 类型

# 爬虫请求头信息
user_agent = 'Mozilla/5.0 (Windows NT 6.1; Win64; x64; rv:48.0) Gecko/20100101 Firefox/48.0'
referer = 'www.genome.jp'
headers = {'User-Agent':user_agent,'Referer':referer}
# 标记该网页中是否有对应的 drug id
str = 'Drug&nbsp;target'
# 为了检测程序，先测试 3 个 target id
# targetID = targetID[0:3]
# print 'targetID'
# print targetID
for id in targetID:
    print('处理 KEGG   数据库（hsa）')
    # 提取 hsa id 对应的所有 drug id
    # DataFrame 格式
    drugID = data.loc[data['Target ID'] == id,'Drug ID']
    id_number = id[3:]       # 提取 hsa 的数字部分
    # hsa id 的网址
    url = 'https://www.genome.jp/dbget-bin/www_bget?hsa:{}'.format(id_number)
    flag, response = get_page(url=url, headers=headers)
    if flag:
        # print('%s 访问网址成功！' %id)
        if str in response.text:
            # 进行内容解析
            re_id = re.compile('href="/dbget-bin/www_bget\?dr:D\d*">(D\d*)</a>')
            items = parse_page(pattern=re_id,response=response.text)
            for did in drugID:
                if did in items:
                    CheckedKEGGhsa.append((did,id))

    else:
        print('%s 访问网址失败！' %id)
        unCheckKEGGhsa.append(id)

# print('unCheckKEGG')
# print unCheckKEGG        # KEGG 读取网页失败的 drug id
# print('CheckedKEGG')
# print CheckedKEGG         # KEGG 已经验证的 (drug id, target id), tuple 类型
# print('unCheckDB')
# print unCheckDB          # DrugBank 网页读取失败的 drug id
# print('CheckedDB')
# print CheckedDB          # DrugBank 已经验证的 (drug id, target id), tuple 类型
# print('unCheckChEMBL')
# print unCheckChEMBL      # ChEMBL 网页读取失败的 drug id
# print('CheckedChEMBL')
# print CheckedChEMBL      # ChEMBL 已经验证的 (drug id, target id), tuple 类型
print('unCheckKEGGhsa')
print unCheckKEGGhsa      # KEGG 网页读取失败的 target id
print('CheckedKEGGhsa')
print CheckedKEGGhsa      # KEGG 已经验证的 (drug id, target id), tuple 类型

df_unkegg = []
df_kegg = []
df_undb = []
df_db = []
df_unchembl = []
df_chembl = []
df_unkegg_hsa = []
df_kegg_hsa = []
# 把非空的结果保存为　.csv
if unCheckKEGG:
    df_unkegg = pd.DataFrame(data=unCheckKEGG,columns=['Drug ID'])
    df_unkegg.to_csv(outName2,index=False)

if CheckedKEGG:
    df_kegg = pd.DataFrame(data=CheckedKEGG,columns=['Drug ID','Target ID'])
    df_kegg['KEGG'] = 'K'          # 增加一列，标记为 K
    print('df kegg')
    print df_kegg
    df_kegg.to_csv(outName1,index=False)

if unCheckDB:
    df_undb = pd.DataFrame(data=unCheckDB,columns=['Drug ID'])
    df_undb.to_csv(outName4,index=False)

if CheckedDB:
    df_db = pd.DataFrame(data=CheckedDB,columns=['Drug ID','Target ID'])
    df_db['DrugBank'] = 'D'         # 增加一列，标记为 D
    print('df db')
    print df_db
    df_db.to_csv(outName3,index=False)

if unCheckChEMBL:
    df_unchembl = pd.DataFrame(data=unCheckChEMBL,columns=['Drug ID'])
    df_unchembl.to_csv(outName6,index=False)

if CheckedChEMBL:
    df_chembl = pd.DataFrame(data=CheckedChEMBL,columns=['Drug ID','Target ID'])
    df_chembl['ChEMBL'] = 'C'
    print('df chembl')
    print df_chembl
    df_chembl.to_csv(outName5,index=False)

if unCheckKEGGhsa:
    df_unkegg_hsa = pd.DataFrame(data=unCheckKEGGhsa,columns=['Target ID'])
    df_unkegg_hsa.to_csv(outName9,index=False)

if CheckedKEGGhsa:
    df_kegg_hsa = pd.DataFrame(data=CheckedKEGGhsa,columns=['Drug ID','Target ID'])
    df_kegg_hsa['KEGGhsa'] = 'K'
    print('df kegg hsa')
    print df_kegg_hsa
    df_kegg_hsa.to_csv(outName8,index=False)

# 考虑到空 list 无法与 DataFrame 进行拼接的，
# 处理方法为将空 list 变成只有 columns=['Drug ID','Target ID'] 的空 DataFrame
# 可以不影响 pd.concat 函数
cols = ['Drug ID','Target ID']
if len(df_chembl) == 0:
    df_chembl = pd.DataFrame([],columns=cols)

if len(df_db) == 0:
    df_db = pd.DataFrame([],columns=cols)

if len(df_kegg) == 0:
    df_kegg = pd.DataFrame([],columns=cols)

if len(df_kegg_hsa) == 0:
    df_kegg_hsa = pd.DataFrame([],columns=cols)

# 3 个 DataFrame 可以进行合并
CheckedInteraction = pd.concat([df_chembl[cols],df_db[cols],df_kegg[cols],df_kegg_hsa[cols]])
# 去除重复
CheckedInteraction = CheckedInteraction.drop_duplicates()
print 'CheckedInteraction'
print CheckedInteraction

# 将 KEGG, DrugBank, ChEMBL数据库所代表的 K, D, C 合并到一个 DataFrame 中
CheckedInteraction = pd.merge(CheckedInteraction,df_chembl,on=cols,how='left')
CheckedInteraction = pd.merge(CheckedInteraction,df_db,on=cols,how='left')
CheckedInteraction = pd.merge(CheckedInteraction,df_kegg,on=cols,how='left')
# 将通过 target ID 在 KEGG 中查找到的作用关系合并到一个 DataFrame 中
CheckedInteraction = pd.merge(CheckedInteraction,df_kegg_hsa,on=cols,how='left')
CheckedInteraction['Checked'] = 'T'
CheckedInteraction.to_csv(outName,index=False)
print 'CheckedInteraction merge'
print CheckedInteraction

# 将查找到的结果合并到原始数据 csv 中
data = pd.merge(data,CheckedInteraction,on=cols,how='left')
data.to_csv(outName7,index=False)
# print outName7
# print data


