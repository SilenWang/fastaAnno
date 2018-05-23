#!/usr/bin/python2
# -*- coding=utf-8 -*-

'''
Author: Silen Wang 
Date: 2018-05-23 22:26:19 
Last Modified by:   Silen Wang 
Last Modified time: 2018-05-23 22:26:19
'''

import sys
import argparse
from argparse import RawTextHelpFormatter
import textwrap


def mark_snp(start, pos, base_list, ref, alt):
    rel_index = pos - start
    # print rel_index
    base = base_list[rel_index]
    if base == ref:  # 使用列表中信息替换原碱基
        base_list[rel_index] = '<b style="color:red">' + alt + '</b>'
        # print base + '<b style="color:red">' + alt + '</b>'
    else:
        print "Warning: base dismatch at %s fas_ref: %s input_ref: %s" % \
            (pos, base, ref)
    return base_list


def mark_indel(start, pos, base_list, ref, alt):
    rel_index = pos - start    
    # print rel_index
    if len(ref) == 1:  # insert 向其中插入多余碱基(会导致行错位)
        base_list[rel_index] = alt[0] + '<b><ins style="color:blue">' + alt[1:] + '</ins></b>'
    else:  # del
        del_count = len(ref) - len(alt)
        base_list[rel_index + 1:rel_index + del_count + 1] =\
            ['<b><del style="color:blue">' + i + '</del></b>'
                for i in base_list[rel_index + 1:rel_index + del_count + 1]]  # 列表生成式挨个替换删除字符
    return base_list


def mark_noDept(start, pos, base_list):
    rel_index = pos - start
    base = base_list[rel_index]
    base_list[rel_index] = '<b style="color: white; background-color: black">' + base + '</b>'
    return base_list


class fas_anno():
    def __init__(self, fas, snp, indel, uncov, out, length, start):
        self.fas = fas
        self.snp = snp
        self.indel = indel
        self.uncov = uncov
        self.out = out
        self.seq_list = []
        self.dis_length = length
        self.start = start

    def mark(self, type):
        if type == "snp":
            with open(self.snp, "r") as snp_file:
                for line in snp_file:
                    if line.startswith("Priority"):
                        head_list = line.strip().split("\t")
                        pos_index = head_list.index("POS")
                        ref_index = head_list.index("REF")
                        alt_index = head_list.index("ALT")
                    else:
                        line_list = line.strip().split("\t")
                        # print line_list
                        self.seq_list = mark_snp(self.start,
                                            int(line_list[pos_index]),
                                            self.seq_list,
                                            line_list[ref_index],
                                            line_list[alt_index])
        elif type == "indel":
            with open(self.indel, "r") as indel_file:
                for line in indel_file:
                    if line.startswith("Priority"):
                        head_list = line.strip().split("\t")
                        pos_index = head_list.index("POS")
                        ref_index = head_list.index("REF")
                        alt_index = head_list.index("ALT")
                    else:
                        line_list = line.strip().split("\t")
                        # print line_list
                        self.seq_list = mark_indel(self.start,
                                            int(line_list[pos_index]),
                                            self.seq_list,
                                            line_list[ref_index],
                                            line_list[alt_index])
        elif type == "uncov":
            # 标注覆盖度为0的位点
            with open(self.uncov, "r") as noDept_file:
                for line in noDept_file:
                    line_list = line.strip().split("\t")
                    self.seq_list = mark_noDept(self.start, int(line_list[1]), self.seq_list)


    def write_seq(self):
        with open(self.out, "w") as out_gene_fas:
            head_str='''
                <!DOCTYPE html>
                <html lang="en">
                <head>
                    <meta charset="UTF-8">
                    <meta name="viewport" content="width=device-width, initial-scale=1.0">
                    <meta http-equiv="X-UA-Compatible" content="ie=edge">
                    <title>Seq_anno_for_mutation</title>

                    <style>
                        body {
                            font-family:monospace;
                        }
                        table tr {
                            line-height: 24px;
                        }
                    </style>
                </head>

                <body>
                <h1>文件说明</h1>
                    <h2>序列表格条目说明:</h2>
                        <p>使用参考基因组: hg19<br>表格分为三列:</p>
                        <table>
                            <tr><td>start_position:</td><td>参考基因组上碱基起始位置</td></tr>
                            <tr><td>sequence_with_anno:</td><td>序列及变异情况</td></tr>
                            <tr><td>end_position:</td><td>参考基因组上碱基终止位置</td></tr>
                        </table>
                    <h2>图例:</h2>
                        <table>
                            <tr><td>点突变(SNP):</td><td>AT<b style="color:red">G</b>C</td></tr>
                            <tr><td>插入(insertion): </td><td>AT<b><ins style="color:blue">GC</ins></b></td></tr>
                            <tr><td>缺失(deletion): </td><td>AT<b><del style="color:blue">GC</del></b></td></tr>
                            <tr><td>测序深度为0X: </td><td>AT<b style="color: white; background-color: black">GC</b></td></tr>
                        </table>

                <hr>
                <table>
                '''
            head_str = textwrap.dedent(head_str)
            out_gene_fas.write(head_str)
            out_gene_fas.write("<tr><th>start_position</th><th>sequence_with_anno</th><th>end_position</th></tr>\n")
            line_count = int(len(self.seq_list) / self.dis_length) + 1  # 算法有隐患, 日后再更改
            # print "line count:" + str(line_count)
            for i in range(line_count):  # 定次循环
                if i != line_count - 1:
                    # out_gene_fas.write("".join(self.seq_list[self.dis_length * (i - 1):self.dis_length * i]) + "<br>\n")  # 写入前self.dis_length个碱基
                    beg_pos = str(self.start + self.dis_length * i)
                    end_pos = str(self.start + self.dis_length * (i+1) - 1)
                    out_gene_fas.write('<tr><td>%s</td><td>%s</td><td>%s</td></tr>\n' %
                        (beg_pos, "".join(self.seq_list[self.dis_length * i:self.dis_length * (i+1) - 1]), end_pos)) # 向表格中写入前self.dis_length个碱基
                else:
                    # print i
                    # print "".join(self.seq_list[self.dis_length * (i - 1):]) + "\n"
                    beg_pos = str(self.start + self.dis_length * i)
                    end_pos = str(self.start + len(self.seq_list) - 1)
                    # out_gene_fas.write("".join(self.seq_list[self.dis_length * (i - 1):]) + "<br>\n")  # 写入最后一行不足self.dis_length个碱基
                    out_gene_fas.write('<tr><td>%s</td><td>%s</td><td>%s</td></tr>\n' %
                        (beg_pos, "".join(self.seq_list[self.dis_length * i:len(self.seq_list) - 1]), end_pos))  # 向表格中写入最后一行不足self.dis_length个碱基
            tail_str="</table>\n</body>\n</html>"
            out_gene_fas.write(tail_str)


    def read_seq(self):
        # 读取完整参考基因组
        with open(self.fas, "r") as ref_gene_fas:
            seq_str = ""
            for seq in ref_gene_fas:
                if not seq.startswith(">"):
                    seq_str += seq.strip()  # 将每行合并成完整序列字符串
            self.seq_list = [i for i in seq_str]  # 列表生成式拆分字符串, 形成参考基因组列表
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description="script to display mutation on refseq", 
            formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        '-r', '--ref',
        help='specify ref genome seq in fasta format',
        type=str,
        required=True
    )
    parser.add_argument(
        '-s', '--snp',
        help='specify snp mutation info tab',
        type=str,
        required=False
    )
    parser.add_argument(
        '-i', '--indel',
        help='specify ref genome seq in fasta format',
        type=str,
        required=False
    )
    parser.add_argument(
        '-o', '--output',
        help='specify the out put file',
        type=str,
        required=True
    )
    parser.add_argument(
        '-u', '--uncov',
        help='specify file that stat uncovered sites',
        type=str,
        required=False
    )
    parser.add_argument(
        '-l', '--length',
        help='specify displayer length in output file',
        type=int,
        # required=True,
        default=150
    )
    parser.add_argument(
        '-S', '--start',
        help='specify pos start for fas file on genome',
        type=int,
        required=True,
        default=1
    )
    args = parser.parse_args()
    fas = fas_anno(
        args.ref, 
        args.snp, 
        args.indel, 
        args.uncov, 
        args.output, 
        args.length, 
        args.start
    )
    fas.read_seq()
    fas.mark("snp")
    fas.mark("indel")
    fas.mark("uncov")
    fas.write_seq()
