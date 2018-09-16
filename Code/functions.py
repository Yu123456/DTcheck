# coding:utf-8
# 2018-06-12

# 函数定义

import requests
import re
import chardet
import numpy as np
import time
import pandas as pd


# 获取网页
def get_page(url,headers):
    # 获取链接 url 的网页内容
    # 返回 (逻辑值，None 或者 content)

    n = 5      # 最多访问网址 5 次
    i = 1
    while i < 5:
        try:
            # 获取网页内容
            r = requests.get(url=url,headers=headers)
            # 通过状态码判断是否获取成功
            if r.status_code == 200:
                r.encoding = chardet.detect(r.content)['encoding']
                return True,r
        except:
            i += 1
        time.sleep(np.random.random() * 5.0)
    return False,None


# 1、Target 栏
def get_target_id(response):
    '''

    :param response:  requests.text 对象
    :return:
    '''

    r_target = re.compile('<th class=.*?"><nobr>Target</nobr></th>')
    result_id = []       # 保存 hsa id
    if re.search(r_target,response):
        # 如果存在 r_target 正则项
        # 先提取 Target 下一行的内容（会包含 Target 中的 HSA）
        r_target_hsa = re.compile('<nobr>Target</nobr></th>\n<td class=.*?"><div(.*?)</div></td>')
        str_hsa = re.findall(r_target_hsa,response)
        # print(str_hsa[0])
        r_hsa = re.compile('<a href=".*?\?hsa:\d+?">\d+?</a>')
        if re.search(r_hsa,str_hsa[0]):
            # 存在 <a href=".*?\?hsa:\d+?">(\d+?)</a> 子串
            r_hsa_id = re.compile('<a href=".*?\?hsa:\d+?">(\d+?)</a>')
            hsa_id = re.findall(r_hsa_id,str_hsa[0])
            # print 'Target'
            # print hsa_id
            for item in hsa_id:
                result_id.append('hsa'+item)

    return result_id

# 2、Pathway 栏
def get_pathway_id(response):
    '''

    :param response:  requests.text 对象
    :return:
    '''
    r_pathway = re.compile('<th class=.*?"><nobr>&nbsp;&nbsp;Pathway</nobr></th>')
    result_id = []  # 保存 hsa id
    if re.search(r_pathway,response):
        # 如果存在 r_pathway 正则项
        # 先提取 Pathway 下一行的内容（会包含 Pathway 中的 hsa）
        r_pathway_hsa = re.compile('"><nobr>&nbsp;&nbsp;Pathway</nobr></th>\n<td class=.*?">'
                                    +'<table style=(.*?)</td></tr></table></td></tr>')
        str_hsa = re.findall(r_pathway_hsa,response)
        r_hsa = re.compile('">hsa\d+?</a>')
        if re.search(r_hsa,str_hsa[0]):
            # 存在 >hsa\d+?</a> 子串
            # 说明有 hsa id
            r_hsa_id = re.compile('">(hsa\d+?)</a>')
            hsa_id = re.findall(r_hsa_id,str_hsa[0])
            # print 'Pathway'
            # print hsa_id
            for item in hsa_id:
                result_id.append(item)
    return result_id


# 3、 Metabolism 栏
def get_metabolism_id(response):
    '''

    :param response:  requests.text 对象
    :return:
    '''
    r_metabolism = re.compile('<th class=.*?"><nobr>Metabolism</nobr></th>')
    result_id = []  # 保存 hsa id
    if re.search(r_metabolism,response):
        # 如果存在 r_metabolism 正则项
        # 先提取 Metabolism 下一行的内容（会包含 Metabolism 中的 hsa）
        r_metabolism_hsa = re.compile('"><nobr>Metabolism</nobr></th>\n<td class=.*?">'
                                    +'<div(.*?)</div></td>',re.S)
        str_hsa = re.findall(r_metabolism_hsa,response)
        r_hsa = re.compile('hsa:\d+?">\d+?</a>')
        if re.search(r_hsa,str_hsa[0]):
            # 存在 hsa:\d+?">\d+?</a> 子串
            # 说明有 hsa id
            r_hsa_id = re.compile('hsa:\d+?">(\d+?)</a>')
            hsa_id = re.findall(r_hsa_id,str_hsa[0])
            # print 'Metabolism'
            # print hsa_id
            for item in hsa_id:
                result_id.append('hsa'+item)
    return result_id

# 4、验证 Interaction 栏
def get_interaction_id(response):
    '''

    :param response: requests.text 对象
    :return:
    '''
    r_interaction = re.compile('<tr><th class=.*?"><nobr>Interaction</nobr></th>')
    result_id = []  # 保存 hsa id
    if re.search(r_interaction,response):
        # 存在 r_interaction 正则项
        # 先提取 Interaction 下一行内容（可能会包含 Interaction 中的 hsa, 必定会存在一张 img,
        # 但是含有 img 行内不会出现 hsa id 信息）
        r_interaction_hsa = re.compile('<tr><th class=.*?"><nobr>Interaction</nobr></th>\n'
                                    +'<td class=(.*?)"><img align=.*?</a></div></td></tr>',re.S)
        str_hsa = re.findall(r_interaction_hsa,response)
        r_hsa = re.compile('hsa:\d+?">\d+?</a>')
        if re.search(r_hsa,str_hsa[0]):
            # 存在 hsa:\d+?">\d+?</a> 子串
            # 说明有 hsa id
            r_hsa_id = re.compile('hsa:\d+?">(\d+?)</a>')
            hsa_id = re.findall(r_hsa_id,str_hsa[0])
            # print 'Interaction'
            # print hsa_id
            for item in hsa_id:
                result_id.append('hsa'+item)
    return result_id

# 5、获取 ChEMBL, DrugBank 中的两个链接
def get_url(response):
    '''

    :param response:  requests.text 对象
    :return:
    '''
    r_url = re.compile('<tr><th class=.*?"><nobr>Other DBs</nobr></th>')
    chembl_url = []
    drugbank_url = []
    if re.search(r_url,response):
        # 存在 r_url 正则项
        # 先提取 Other DBs 下一行内容
        r_url_dbs = re.compile('<tr><th class=.*?"><nobr>Other DBs</nobr></th>\n'
                               +'<td class=(.*?</a>)</td></tr></table></td></tr>\n')
        str_url = re.findall(r_url_dbs,response)

        # 处理 ChEMBL 链接
        r_chembl = re.compile('<td valign="top"><nobr>ChEMBL:&nbsp;</nobr>'
                              +'</td>(.*?</a>)</td></tr></table>')
        str_chembl = re.findall(r_chembl,str_url[0])
        if str_chembl:
            # str_chembl 非空，说明有 ChEMBL 链接
            r_url_chembl = re.compile('<a href="(http.*?/CHEMBL\d*?)"'
                                      +'>CHEMBL\d*?</a>')
            urls = re.findall(r_url_chembl,str_chembl[0])
            chembl_url.extend(urls)

        # 处理 DrugBank 链接
        r_drugbank = re.compile('<td valign="top"><nobr>DrugBank:&nbsp;</nobr></td>'
                                +'<td>(.*?</a>)</td></tr></table>')
        str_drugbank = re.findall(r_drugbank,str_url[0])
        if str_drugbank:
            # str_drugbank 非空，说明有 DrugBank 链接
            r_url_drugbank = re.compile('<a href="(http.*?/DB\d+?)">DB\d+?</a>')
            urls = re.findall(r_url_drugbank,str_drugbank[0])
            if urls:
                drugbank_url.extend(urls)

    return chembl_url,drugbank_url


# DrugBank 数据库中获取该 drug 所能对应的
# Targets , Enzymes, Carriers, Transporters 的 Uniprot ID
# 尽管有些 drug 不会对应所有 4 类，因为只需要 Uniprot ID，可以从 Targets 开始
# 截获网页内容，提取 Uniprot ID 即可
# 注意，默认所有 drug 都会有 Targets 项，如果运行中出错，再做修改
def get_uniprot_id(response):
    '''

    :param response:  requests.text 对象
    :return:
    '''
    # 1、Targets 栏
    # 2、Enzymes 栏
    # 3、Carriers 栏
    # 4、Transporters 栏

    uniprot_id = []     # 返回的 uniprot id
    # 4 个栏目下的 uniprot id 的模式是一致的
    r_id = re.compile('<a target="_blank" href="http\w*?://www.uniprot.org/'
                      + 'uniprot/\w+?\d+?">(\w+?\d+?)</a>')

    r_targets = re.compile('<h3 id="targets".*?>Targets</h3>')
    if re.search(r_targets,response):
        # 存在 Targets 栏
        r_uniprot_id = re.compile('<div class="bond-list-container targets">'
                                  + '<h3 id="targets".*?>Targets</h3>'
                                  + '<div.*?<p class="bt-2" id="drug-meta">'
                                  + 'Drug created on .*? </p></div>')
        str_uniprot_id = re.search(r_uniprot_id,response)
        if str_uniprot_id:
            u_id = re.findall(r_id,str_uniprot_id.group())
            if u_id:
                uniprot_id.extend(u_id)
        return uniprot_id

    r_enzymes = re.compile('<h3 id="targets".*?>Enzymes</h3>')
    if re.search(r_enzymes,response):
        # 存在 Enzymes 栏
        r_uniprot_id = re.compile('<div class="bond-list-container enzymes">'
                                  + '<h3 id="enzymes".*?>Enzymes</h3>'
                                  + '<div.*?<p class="bt-2" id="drug-meta">'
                                  + 'Drug created on .*? </p></div>')
        str_uniprot_id = re.search(r_uniprot_id, response)
        if str_uniprot_id:
            u_id = re.findall(r_id, str_uniprot_id.group())
            if u_id:
                uniprot_id.extend(u_id)
        return uniprot_id

    r_carriers = re.compile('<h3 id="carriers".*?>Carriers</h3>')
    if re.search(r_carriers,response):
        # 存在 Enzymes 栏
        r_uniprot_id = re.compile('<div class="bond-list-container carriers">'
                                  + '<h3 id="carriers".*?>Carriers</h3>'
                                  + '<div.*?<p class="bt-2" id="drug-meta">'
                                  + 'Drug created on .*? </p></div>')
        str_uniprot_id = re.search(r_uniprot_id, response)
        if str_uniprot_id:
            u_id = re.findall(r_id, str_uniprot_id.group())
            if u_id:
                uniprot_id.extend(u_id)
        return uniprot_id

    r_transporters = re.compile('<h3 id="transporters".*?>Transporters</h3>')
    if re.search(r_carriers, response):
        # 存在 Enzymes 栏
        r_uniprot_id = re.compile('<div class="bond-list-container transporters">'
                                  + '<h3 id="transporters".*?>Transporters</h3>'
                                  + '<div.*?<p class="bt-2" id="drug-meta">'
                                  + 'Drug created on .*? </p></div>')
        str_uniprot_id = re.search(r_uniprot_id, response)
        if str_uniprot_id:
            u_id = re.findall(r_id, str_uniprot_id.group())
            if u_id:
                uniprot_id.extend(u_id)
        return uniprot_id

    # 如果所有都不符合，直接返回空 list
    return uniprot_id


# 把 uniprot id 转换成 KEGG id (hsa id)
def IDmapping(u_id):
    '''
    1、输入要保证是无重复 uniprot id, 重复的只返回一个
       ID mapping 每一页最多显示 25 个 id，超过该数就涉及到翻页
    2、uploadQuery: 值需要合并成一个字符串，以分号分隔，
    :param u_id: uniprot id, 传入的是 list
    :return: 返回 KEGG id
    '''
    hsa_id = []      # 返回 hsa id

    max_n = 25
    ns = len(u_id)     # 需要查询 id 的个数
    if ns > max_n:
        n_splits = ns // max_n    # 至少转换次数
        fold_sizes = max_n * np.ones(n_splits,dtype=np.int)
        if ns - max_n * n_splits > 0:
            fold_sizes = np.append(fold_sizes,ns - max_n * n_splits)

        current = 0
        for fold_size in fold_sizes:
            start,stop = current, current+fold_size
            batch_id = u_id[start:stop]
            batch_id = ";".join(batch_id)
            result_id = get_mapping(batch_id)
            if result_id:
                hsa_id.extend(result_id)
    else:
        batch_id = ";".join(u_id)
        result_id = get_mapping(batch_id)
        if result_id:
            hsa_id.extend(result_id)

    return hsa_id





# 提交 FormTable 获取转换后的 KEGG id
def get_mapping(str_id):
    '''

    :param str_id:
    :return:
    '''

    url = 'https://www.uniprot.org/uploadlists/'
    # 爬虫请求头信息
    user_agent = 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:49.0) Gecko/20100101 Firefox/49.0'
    referer = 'https://www.uniprot.org/uploadlists/'
    headers = {'User-Agent': user_agent, 'Referer': referer}
    formtable = {'landingPage': 'false', 'uploadQuery': str_id, 'from': 'ACC,ID', 'to': 'KEGG_ID'}

    n = 5  # 最多访问网址 5 次
    i = 1
    hsa_id = []     # 返回的 hsa id
    while i < 5:
        try:
            # 获取网页内容
            response = requests.post(url=url, data=formtable,headers=headers)
            # 通过状态码判断是否获取成功
            if response.status_code == 200:
                response.encoding = chardet.detect(response.content)['encoding']
                # 抓取网页中转换后的 hsa
                r_FromTo = re.compile('<thead><tr><th>From</th><th>To</th></tr></thead>'
                                      +'<tbody>.*?</td></tr></tbody></table>')
                str_hsa = re.search(r_FromTo,response.text)
                if str_hsa:
                    # str_hsa 非空
                    r_hsa = re.compile('<a href="http.*?hsa:\d+?">hsa:(\d+?)</a>')
                    id = re.findall(r_hsa,str_hsa.group())
                    for item in id:
                        hsa_id.append('hsa'+item)

                return hsa_id
        except:
            i += 1
    print('转换失败！')

    return hsa_id


# 从 ChEMBL 中获得 uniprot id, 然后转成 KEGG ID
def get_hsaid_chembl(url):
    '''

    :param url:
    :return:
    '''

    hsa_id = []               # 返回转换后的 id
    uniprot_id = []           # 查找到的 uniprot id

    # 获取 url 网页
    # 爬虫请求头信息
    user_agent = 'Mozilla/5.0 (Windows NT 6.1; WOW64; rv:49.0) Gecko/20100101 Firefox/49.0'
    referer = 'www.ebi.ac.uk'
    headers = {'User-Agent': user_agent, 'Referer': referer}
    flag, response = get_page(url=url,headers=headers)

    if flag:
        # 获取网页成功
        r_url = re.compile('<h2 style=.*?>Mechanism of Action'
                       + '</h2>\n<table.*?</tr>\n</table>', re.S)
        if re.search(r_url,response.text):
            # 存在 r_url 正则项
            # 实际上是 https://www.ebi.ac.uk/chembldb/compound/inspect/CHEMBL***， compound 类网址
            new_url = get_fromMechanismofAction(response=response.text)
            for item in new_url:
                newflag,newresponse = get_page(url=item,headers=headers)
                if newflag:
                    # 获取网页成功（含 Target Components ）
                    cur_u_id = get_fromTargetComponents(response=newresponse.text)
                    if cur_u_id:
                        uniprot_id.extend(cur_u_id)

        # 查询是否有 Target Components
        # 实际上是 https://www.ebi.ac.uk/chembl/target/inspect/CHEMBL***, target 类网址
        r_u_id = re.compile('<h2 style=.*?>Target Components'
                            + '</h2>\n<table.*?\n</tr>\n</table>', re.S)
        str_u_id = re.search(r_u_id,response.text)
        if str_u_id:
            cur_u_id = get_fromTargetComponents(response=response.text)
            if cur_u_id:
                uniprot_id.extend(cur_u_id)

    # 如果 uniprot_id 非空，就转换成  KEGG ID
    if uniprot_id:
        uniprot_id = np.unique(uniprot_id)
        result_u_id = IDmapping(u_id=uniprot_id)
        if result_u_id:
            hsa_id.extend(result_u_id)

    return hsa_id




# 从 Mechanism of Action 提取 ChEMBL Target 链接
def get_fromMechanismofAction(response):
    '''

    :param response:  requests.text 对象
    :return:
    '''

    target_url = []      # 返回 ChEMBL Target 链接
    url_start = 'https://www.ebi.ac.uk'
    r_url = re.compile('<h2 style=.*?>Mechanism of Action'
                       +'</h2>\n<table.*?</tr>\n</table>',re.S)
    str_url = re.search(r_url,response)
    if str_url:
        # 如果存在 url
        r_chembl_url = re.compile('<a href=\'(/chembl/target/inspect/CHEMBL\d+?)\'>')
        str_chembl_url = re.findall(r_chembl_url,str_url.group())
        if str_chembl_url:
            for item in str_chembl_url:
                target_url.append(url_start+item)
    return target_url


# 从 Target Components 提取 uniprot id
def get_fromTargetComponents(response):
    '''

    :param response:  requests.text 对象
    :return:
    '''

    uniprot_id = []        # 返回 uniprot id
    r_u_id = re.compile('<h2 style=.*?>Target Components'
                        +'</h2>\n<table.*?\n</tr>\n</table>',re.S)
    str_u_id = re.search(r_u_id,response)
    # print 'search(r_u_id)'
    # print str_u_id
    if str_u_id:
        r_u = re.compile('<a href="http://www.uniprot.org/uniprot'
                         +'/\w*?\d+?">(\w*?\d+?)</a>')
        str_u = re.findall(r_u,str_u_id.group())
        if str_u:
            uniprot_id.extend(str_u)

    return uniprot_id

# 核查药靶作用关系，返回 tuple 作用对构成的 list
# 如果没有作用关系对，返回空 list
def check_interaction(drugid,hsa_predict,hsa_db):
    '''
    核查 hsa_predict 中的 id 是否在 hsa_db 中
    :param drugid: drug id
    :param hsa_predict:   预测与 drug id 有作用关系的 hsa id
    :param hsa_db:  数据库中查询到与 drug id 有作用关系的 hsa id
    :return result: list 类型，包含 tuple 对
    '''
    # print 'check_interaction'
    result = []
    # 转换成 set ，然后求交集
    # 只有 hsa_predict 与 hsa_db 均非空时，才进行交集运算，否则返回一个空 list
    if len(hsa_db) and len(hsa_predict):
        s1 = set(hsa_predict)
        s2 = set(hsa_db)
        s = s1.intersection(s2)
        for item in s:
            result.append((drugid,item))

    return result


# 合并两个 DataFrame
def df_merge(df_left,df_right,on,how):
    '''

    :param df_left:  需要合并的数据框
    :param df_right:  需要合并的数据框
    :param on:   连接的列名
    :param how:  连接方式(只处理 left , right 两种方式)
    :return:
    '''

    df = []
    if len(df_left) and len(df_right):
        # 两个均非空
        df = pd.merge(df_left,df_right,on=on,how=how)
    elif len(df_left):
        # df_left 非空
        # 注意，进入该分支表明 df_left 非空而 df_right 一定为空
        # 并且此时第一个 if 条件为假，即必有一个为空，或者两者均为空
        if on == 'left':
            # 左连接方式
            df = df_left
        else:
            # 右连接方式
            df = df_right
    else:
        # 此时 df_right 可能为空，也可能非空
        # 但是 df_left 一定为空
        if on =='left':
            # 左连接方式
            df = []
        else:
            # 右连接方式
            df = df_right

    return df


# 拼接两个 DataFrame
def df_concat(df1,df2,col_name):
    '''
    只对 df1, df2 按照列名拼接
    :param df1:  DataFrame
    :param df2:   DataFrame
    :param col_name:  ['name1','name2'] 拼接的列名，也只返回包含该列名的 DataFrame
    :return:
    '''

    df = []
    if len(df1) and len(df2):
        # 两个均非空
        df = pd.concat(df1[col_name],df2[col_name])
    elif len(df1):
        # df1 非空，
        # 此时说明 df2 一定为空
        df = df1[col_name]
    elif len(df2):
        # df2 非空
        # 此时说明 df1 一定为空
        df = df2[col_name]

    if len(df):
        # 如果 df 非空，对其进行去重复处理
        print('去除重复')
        df = df.drop_duplicates()

    return df

# hsa Check 时的函数
# 解析网页，返回需要内容
def parse_page(pattern,response):
    '''

    :param pattern: 需要匹配的模式
    :param response: requests.text 格式
    :return:
    '''

    return re.findall(pattern,response)