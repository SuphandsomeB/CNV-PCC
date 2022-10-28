import numpy as np
import pysam
import os
import sys
import pandas as pd
from numba import njit
import subprocess
from pyod.models.pca import PCA
from imblearn.under_sampling import RandomUnderSampler

def get_chrlist(filename):

    samfile = pysam.AlignmentFile(filename, "rb", ignore_truncation=True)
    List = samfile.references
    chrList = []
    for i in range(len(List)):
        chr = str(List[i]).strip('chr')
        if chr.isdigit():
            chrList.append(int(chr))
    return chrList

def get_RC(filename ,ReadCount, Mapq):

    qname, flag, q_mapq, cigar = [], [], [], []
    pos, direction, isize, qlen = [], [], [], []
    samfile = pysam.AlignmentFile(filename, "rb", ignore_truncation=True)

    for line in samfile.fetch(chr_name):
        if line.reference_name:
            chr = int(line.reference_name.strip('chr'))
            if chr == chr_num:
                posList = line.positions
                ReadCount[posList] += 1
                Mapq[posList] += line.mapq
                if line.cigarstring != None:
                    cigarstring = str(line.query_length) + 'M'
                    if line.cigarstring != cigarstring and line.mapq >= 10:
                        qname.append(line.qname)
                        flag.append(line.flag)
                        q_mapq.append(line.mapq)
                        cigar.append(line.cigarstring)
                        pos.append(line.pos)
                        direction.append(line.is_read1)
                        isize.append(line.isize)
                        qlen.append(line.qlen)
    SR = [*zip(qname, flag, direction, q_mapq, cigar, pos, qlen)]
    SR = pd.DataFrame(SR, columns=['name', 'flag', 'dir', 'mapq', 'cigar', 'pos', 'len'])
    Bpoint = get_breakpoint2(SR)

    return ReadCount, Bpoint, Mapq

def get_breakpoint(read_cigar,read_pos,read_len):
    breakpoint = []
    for i in range(len(read_cigar)):
        if 'S' in read_cigar[i] and 'M' in read_cigar[i]:
            S_index = read_cigar[i].index('S')
            M_index = read_cigar[i].index('M')
            if S_index < M_index:
                breakpoint.append(read_pos[i])
            elif S_index > M_index:
                breakpoint.append(read_pos[i] + read_len[i])
        elif 'H' in read_cigar[i] and 'M' in read_cigar[i]:
            H_index = read_cigar[i].index('H')
            M_index = read_cigar[i].index('M')
            if H_index < M_index:
                breakpoint.append(read_pos[i])
            elif H_index > M_index:
                breakpoint.append(read_pos[i] + read_len[i])
        else:
            continue
    return breakpoint

def get_breakpoint2(SR):
    breakpoint1 = []
    breakpoint2 = [0,chrLen[chr_num-1]-1]
    for name, group in SR.groupby('name'):
        num_g = group['name'].count()
        if num_g == 1:
            read_cigar = np.array(group['cigar'])
            read_pos = np.array(group['pos'])
            read_len = np.array(group['len'])
            breakpoint1 += get_breakpoint(read_cigar, read_pos, read_len)
        elif num_g >= 2:
            index1 = np.array(group['dir']) == True
            read1_count = np.sum(index1)
            index2 = np.array(group['dir']) == False
            read2_count = np.sum(index2)

            if read1_count > 1:
                read1_cigar = np.array(group['cigar'][index1])
                read1_pos = np.array(group['pos'][index1])
                read1_len = np.array(group['len'][index1])
                breakpoint2 += get_breakpoint(read1_cigar, read1_pos, read1_len)

            if read2_count > 1:
                read2_cigar = np.array(group['cigar'][index2])
                read2_pos = np.array(group['pos'][index2])
                read2_len = np.array(group['len'][index2])
                breakpoint2 += get_breakpoint(read2_cigar, read2_pos, read2_len)

    breakpoint2 = sorted(list(set(breakpoint2)))
    return breakpoint2

def read_ref_file(filename, ref):
    # read  reference.file
    # input：.fasta file
    if os.path.exists(filename):
        print("Read reference file: " + str(filename))
        with open(filename, 'r') as f:
            line = f.readline()
            chr_name = line.strip('>').strip()
            chr_num = int(line.strip('>chr'))
            for line in f:
                linestr = line.strip()
                ref[chr_num - 1] += linestr
    else:
        print("Warning: can not open " + str(filename) + '\n')
    return ref,chr_num,chr_name

def ReadDepth(ReadCount, ref, pos):

    start = pos[:len(pos) - 1]
    end = pos[1:]
    for i in range(len(pos) - 1):
        if end[i] - start[i] < 500:
            pos.remove(end[i])
    pos = np.array(pos)
    start = pos[:len(pos) - 1]
    end = pos[1:]
    length = end - start
    with open('pos.txt', 'w') as f:
        for i in range(len(start)):
            linestrlist = ['1', '1', str(start[i]), str(end[i]-1), str(length[i])]
            f.write('\t'.join(linestrlist) + '\n')
    bin_start,bin_end,bin_len,s = re_segfile('pos.txt', 'bin.txt', binSize)

    binNum = len(bin_start)
    bin_RD = np.full(binNum, 0.0)
    bin_GC = np.full(binNum, 0)
    bin_gc = np.full(binNum,0.0)
    for i in range(binNum):
        bin_RD[i] = np.mean(ReadCount[bin_start[i]:bin_end[i]])
        cur_ref = ref[bin_start[i]:bin_end[i]]
        N_count = cur_ref.count('N') + cur_ref.count('n')
        if N_count == 0:
            gc_count = cur_ref.count('C') + cur_ref.count('c') + cur_ref.count('G') + cur_ref.count('g')
        else:
            bin_RD[i] = -10000
            gc_count = 0
        bin_GC[i] = int(round(gc_count / bin_len[i], 3) * 1000)
        bin_gc[i] = round(gc_count / bin_len[i],2)
    bin_end -= 1

    index = bin_RD > 0
    bin_RD = bin_RD[index]
    bin_GC = bin_GC[index]
    bin_gc = bin_gc[index]
    bin_len = bin_len[index]
    bin_start = bin_start[index]
    bin_end = bin_end[index]
    bin_RD = gc_correct(bin_RD, bin_GC)
    return bin_start,bin_end,bin_len, bin_RD, bin_gc


def gc_correct(RD, GC):

    bincount = np.bincount(GC)
    global_rd_ave = np.mean(RD)
    for i in range(len(RD)):
        if bincount[GC[i]] < 2:
            continue
        mean_RD = np.mean(RD[GC == GC[i]])
        RD[i] = global_rd_ave * RD[i] / mean_RD
    return RD


def prox_tv1d(step_size: float, w: np.ndarray) -> np.ndarray:

    if w.dtype not in (np.float32, np.float64):
        raise ValueError('argument w must be array of floats')
    w = w.copy()
    output = np.empty_like(w)
    _prox_tv1d(step_size, w, output)
    return output

@njit
def _prox_tv1d(step_size, input, output):
    """low level function call, no checks are performed"""
    width = input.size + 1
    index_low = np.zeros(width, dtype=np.int32)
    slope_low = np.zeros(width, dtype=input.dtype)
    index_up  = np.zeros(width, dtype=np.int32)
    slope_up  = np.zeros(width, dtype=input.dtype)
    index     = np.zeros(width, dtype=np.int32)
    z         = np.zeros(width, dtype=input.dtype)
    y_low     = np.empty(width, dtype=input.dtype)
    y_up      = np.empty(width, dtype=input.dtype)
    s_low, c_low, s_up, c_up, c = 0, 0, 0, 0, 0
    y_low[0] = y_up[0] = 0
    y_low[1] = input[0] - step_size
    y_up[1] = input[0] + step_size
    incr = 1

    for i in range(2, width):
        y_low[i] = y_low[i-1] + input[(i - 1) * incr]
        y_up[i] = y_up[i-1] + input[(i - 1) * incr]

    y_low[width-1] += step_size
    y_up[width-1] -= step_size
    slope_low[0] = np.inf
    slope_up[0] = -np.inf
    z[0] = y_low[0]

    for i in range(1, width):
        c_low += 1
        c_up += 1
        index_low[c_low] = index_up[c_up] = i
        slope_low[c_low] = y_low[i]-y_low[i-1]
        while (c_low > s_low+1) and (slope_low[max(s_low, c_low-1)] <= slope_low[c_low]):
            c_low -= 1
            index_low[c_low] = i
            if c_low > s_low+1:
                slope_low[c_low] = (y_low[i]-y_low[index_low[c_low-1]]) / (i-index_low[c_low-1])
            else:
                slope_low[c_low] = (y_low[i]-z[c]) / (i-index[c])

        slope_up[c_up] = y_up[i]-y_up[i-1]
        while (c_up > s_up+1) and (slope_up[max(c_up-1, s_up)] >= slope_up[c_up]):
            c_up -= 1
            index_up[c_up] = i
            if c_up > s_up + 1:
                slope_up[c_up] = (y_up[i]-y_up[index_up[c_up-1]]) / (i-index_up[c_up-1])
            else:
                slope_up[c_up] = (y_up[i]-z[c]) / (i-index[c])

        while (c_low == s_low+1) and (c_up > s_up+1) and (slope_low[c_low] >= slope_up[s_up+1]):
            c += 1
            s_up += 1
            index[c] = index_up[s_up]
            z[c] = y_up[index[c]]
            index_low[s_low] = index[c]
            slope_low[c_low] = (y_low[i]-z[c]) / (i-index[c])
        while (c_up == s_up+1) and (c_low>s_low+1) and (slope_up[c_up]<=slope_low[s_low+1]):
            c += 1
            s_low += 1
            index[c] = index_low[s_low]
            z[c] = y_low[index[c]]
            index_up[s_up] = index[c]
            slope_up[c_up] = (y_up[i]-z[c]) / (i-index[c])

    for i in range(1, c_low - s_low + 1):
        index[c+i] = index_low[s_low+i]
        z[c+i] = y_low[index[c+i]]
    c = c + c_low-s_low
    j, i = 0, 1
    while i <= c:
        a = (z[i]-z[i-1]) / (index[i]-index[i-1])
        while j < index[i]:
            output[j * incr] = a
            output[j * incr] = a
            j += 1
        i += 1
    return


@njit
def prox_tv1d_cols(stepsize, a, n_rows, n_cols):
    """apply prox_tv1d along columns of the matri a
    """
    A = a.reshape((n_rows, n_cols))
    out = np.empty_like(A)
    for i in range(n_cols):
        _prox_tv1d(stepsize, A[:, i], out[:, i])
    return out.ravel()


@njit
def prox_tv1d_rows(stepsize, a, n_rows, n_cols):
    """apply prox_tv1d along rows of the matri a
    """
    A = a.reshape((n_rows, n_cols))
    out = np.empty_like(A)
    for i in range(n_rows):
        _prox_tv1d(stepsize, A[i, :], out[i, :])
    return out.ravel()

def Read_seg_file(binstart,binlen,binend,bingc):
    """
    read segment file (Generated by DNAcopy.segment)
    seg file: col, chr, start, end, num_mark, seg_mean
    """
    seg_start = []
    seg_end = []
    seg_len = []
    location = []
    seg_pos = []
    for i in range(len(binstart) - 1):
        if binstart[i] + binlen[i] != binstart[i + 1]:
            location.append(i+1)
    count = 0
    with open("seg", 'r') as f1,\
            open('seg2.txt', 'w') as f2:
        for line in f1:
            linestrlist = line.strip().split('\t')
            start = int(linestrlist[2])
            end = int(linestrlist[3])
            seg_pos.append([])
            seg_pos[-1].append(start)
            seg_pos[-1].append(end)
            for j in location:
                if j >= seg_pos[-1][0] and j < seg_pos[-1][-1]:
                    seg_pos[-1].append(j)
                    seg_pos[-1].append(j + 1)
                    seg_pos[-1] = sorted(seg_pos[-1])
            for k in range(0,len(seg_pos[count]),2):
                start = seg_pos[count][k] - 1
                end = seg_pos[count][k+1]-1
                linestrlist[2] = str(binstart[seg_pos[count][k] - 1])
                seg_start.append(seg_pos[count][k]-1)
                linestrlist[3] = str(binend[seg_pos[count][k+1]-1])#
                linestrlist[4] = str(np.sum(binlen[seg_pos[count][k] - 1:seg_pos[count][k+1]]))
                seg_end.append(seg_pos[count][k+1]-1)
                seg_len.append(np.sum(binlen[start:end + 1]))
                linestrlist.append('')
                linestrlist[6] = (str(np.mean(bingc[start:end+1])))
                f2.write('\t'.join(linestrlist) + '\n')
            count += 1
    reseg_Start,reseg_End,reseg_Len,reseg_gc = re_segfile('seg2.txt','reseg.txt',reseg_len)
    return reseg_Start, reseg_End, reseg_Len, reseg_gc

def PCC(data):
    rdmq_1 = np.array(data[['rd','gc','mq']])
    clf = PCA(n_components=1)
    clf.fit(rdmq_1)
    scores_1 = clf.decision_scores_
    data['scores'] = scores_1
    threshold_1 = Otsu(scores_1)
    CNV_1 = get_CNV(data,threshold_1)
    CNV_stage1 = combineCNV(CNV_1)
    new_data = data.drop(CNV_1.index)
    new_data.index = range(new_data.shape[0])
    rdmq_2 = np.array(new_data[['rd','gc','mq']])
    clf.fit(rdmq_2)
    scores_2 = clf.decision_scores_
    new_data['scores'] = scores_2
    threshold_2 = Otsu(scores_2)
    CNV_2 = get_CNV(new_data, threshold_2)
    CNV_stage2 = combineCNV(CNV_2)
    allCNV = pd.concat([CNV_1, CNV_2]).sort_values(by='start').reset_index(drop=True)
    final_CNV = combineCNV(allCNV)

    return final_CNV

def resegment_RD(RD, MQ, start, end):

    reseg_RD = np.full(len(start),0.0)
    reseg_MQ = np.full(len(start),0.0)
    reseg_start = []
    reseg_end = []

    for i in range(len(start)):
        reseg_RD[i] = np.mean(RD[start[i]:end[i]])
        reseg_MQ[i] = np.mean(MQ[start[i]:end[i]])
        reseg_start.append(start[i] + 1)
        reseg_end.append(end[i])

    return reseg_RD, reseg_MQ, reseg_start, reseg_end

def re_segfile(filname,savefile,reseg_length):
    with open(filname,'r') as f1, \
            open(savefile,'w') as f2:
        for line in f1:
            linestrlist = line.strip().split('\t')
            length = int(linestrlist[4])
            if length > reseg_length:
                l = length % reseg_length
                if l:
                    if l >= binSize:
                        reseg_num = length // reseg_length + 1
                        for i in range(reseg_num):
                            if i == 0:
                                linestrlist[3] = str(int(linestrlist[2]) + reseg_length - 1)
                                linestrlist[4] = str(reseg_length)
                                f2.write('\t'.join(linestrlist) + '\n')
                            elif i+1 != reseg_num:
                                linestrlist[2] = str(int(linestrlist[2]) + 1 * reseg_length)
                                linestrlist[3] = str(int(linestrlist[2]) + reseg_length - 1 )
                                linestrlist[4] = str(reseg_length)
                                f2.write('\t'.join(linestrlist) + '\n')
                            else:
                                linestrlist[2] = str(int(linestrlist[2]) + 1 * reseg_length)
                                linestrlist[3] = str(int(linestrlist[2])  + l - 1)
                                linestrlist[4] = str(l)
                                f2.write('\t'.join(linestrlist) + '\n')

                    else:
                        bin_num = length // reseg_length
                        for i in range(bin_num):
                            if i == 0:
                                if i+1 != bin_num:
                                    linestrlist[3] = str(int(linestrlist[2]) + reseg_length - 1)
                                    linestrlist[4] = str(reseg_length)
                                else:
                                    linestrlist[3] = str(int(linestrlist[2]) + reseg_length - 1 + l)
                                    linestrlist[4] = str(reseg_length + l)
                                f2.write('\t'.join(linestrlist) + '\n')
                            elif i+1 != bin_num:
                                linestrlist[2] = str(int(linestrlist[2]) + 1 * reseg_length)
                                linestrlist[3] = str(int(linestrlist[2]) + reseg_length - 1 )
                                linestrlist[4] = str(reseg_length)
                                f2.write('\t'.join(linestrlist) + '\n')
                            else:
                                linestrlist[2] = str(int(linestrlist[2]) + 1 * reseg_length)
                                linestrlist[3] = str(int(linestrlist[2]) + reseg_length - 1 + l)
                                linestrlist[4] = str(reseg_length + l)
                                f2.write('\t'.join(linestrlist) + '\n')

                else:
                    bin_num = length // reseg_length
                    for i in range(bin_num):
                        if i==0:
                            linestrlist[3] = str(int(linestrlist[2]) + reseg_length - 1)
                            linestrlist[4] = str(reseg_length)
                            f2.write('\t'.join(linestrlist) + '\n')
                        else:
                            linestrlist[2] = str(int(linestrlist[2]) + 1 * reseg_length)
                            linestrlist[3] = str(int(linestrlist[2]) + reseg_length - 1)
                            linestrlist[4] = str(reseg_length)
                            f2.write('\t'.join(linestrlist) + '\n')
            else:
                f2.write('\t'.join(linestrlist) + '\n')

    tran_start = []
    tran_end = []
    tran_len = []
    tran_gc = []
    with open(savefile, 'r') as f:
        for line in f:
            linestrinfo = line.strip().split('\t')
            tran_start.append(int(linestrinfo[2]))
            tran_end.append(int(linestrinfo[3]) + 1)
            tran_len.append(int(linestrinfo[4]))
            if len(linestrinfo) > 5:
                tran_gc.append(round(float(linestrinfo[6]),2))

    tran_start = np.array(tran_start)
    tran_end = np.array(tran_end)
    tran_len = np.array(tran_len)
    os.remove(filname)
    os.remove(savefile)
    return tran_start, tran_end, tran_len,tran_gc


def get_newbins(new_data):
    new_chr = np.array(new_data['chr'])
    new_start = np.array(new_data['start'])
    new_end = np.array(new_data['end'])
    new_rd = np.array(new_data['rd'])
    new_gc = np.array(new_data['gc'])
    new_mq = np.array(new_data['mq'])

    return new_chr,new_start,new_end,new_rd,new_gc,new_mq

def Otsu(S):
    S = np.round(S,2)
    min_S = np.min(S)
    median_S = np.median(S)
    lower_S = np.quantile(S,0.35,interpolation='lower')
    higer_S = np.quantile(S,0.85,interpolation='higher')
    if(lower_S == min_S):
        lower_S += 0.1
    final_threshold = median_S
    max_var = 0.0
    D_labels = np.full(len(S), 0)
    for i in np.arange(lower_S,higer_S,0.01):
        cur_threshold = round(i,2)
        D0_index = (S < cur_threshold)
        D1_index = (S >= cur_threshold)

        D_labels[D0_index] = 0
        D_labels[D1_index] = 1
        D0 = S[D0_index]
        D1 = S[D1_index]
        S_resample = S.reshape(-1,1)

        new_D,new_label = RandomUnderSampler(random_state=42).fit_resample(S_resample,D_labels)
        new_D0 = new_D.ravel()[new_label==0]
        new_D1 = new_D.ravel()[new_label==1]

        D0_mean = np.mean(new_D0)
        D1_mean = np.mean(new_D1)
        p0 = len(D0)/(len(D0)+len(D1))
        p1 = (1 - p0)
        S_mean = p0*D0_mean + p1*D1_mean
        cur_var = p0*(D0_mean - S_mean)**2 + p1*(D1_mean - S_mean)**2
        if cur_var > max_var:
            final_threshold = cur_threshold
            max_var = cur_var

    return final_threshold

def get_CNV(data,threshold):

    CNVindex = data[np.round(data['scores'],2) >= threshold].index
    Normalindex = data[np.round(data['scores'],2) < threshold].index
    Normalmean = (data['rd'].iloc[Normalindex]).mean()
    gc_mean = (data['gc'].iloc[Normalindex]).mean()
    base = Normalmean * 0.25
    base2 = Normalmean * 0.5
    r1 = data['rd'].iloc[CNVindex] < Normalmean - base2
    r2 = data['rd'].iloc[CNVindex] > Normalmean + base
    real_CNVindex_loss = CNVindex[r1.values]
    real_CNVindex_gain = CNVindex[r2.values]
    cnv_loss = data.iloc[real_CNVindex_loss]
    cnv_gain = data.iloc[real_CNVindex_gain]
    CNVtype_loss = np.full(cnv_loss.shape[0], 'loss')
    cnv_loss.insert(6, 'type', CNVtype_loss)
    CNVtype_gain = np.full(cnv_gain.shape[0], 'gain')
    cnv_gain.insert(6, 'type', CNVtype_gain)
    allCNV = pd.concat([cnv_loss, cnv_gain])
    CNV_length = allCNV['end'] - allCNV['start'] + 1
    allCNV.insert(3, 'length', CNV_length)

    return allCNV

def combineCNV(CNVdata):

    CNVchr,CNVstart,CNVend,CNVRD,CNVgc,CNVmq = get_newbins(CNVdata)
    CNVlen = CNVend - CNVstart + 1
    typeCNV = np.array(CNVdata['type'])
    for i in range(len(CNVRD) - 1):
        if typeCNV[i] == typeCNV[i + 1]:
            len_n = CNVstart[i + 1] - CNVend[i] - 1
            if len_n / (CNVlen[i] + CNVlen[i + 1] + len_n) == 0:
                CNVstart[i + 1] = CNVstart[i]
                CNVlen[i + 1] = CNVend[i + 1] - CNVstart[i + 1] + 1
                typeCNV[i] = 0
                CNVRD[i + 1] = (CNVRD[i] + CNVRD[i+1]) / 2
                CNVmq[i + 1] = (CNVmq[i] + CNVmq[i + 1]) / 2
                CNVgc[i + 1] = (CNVgc[i] + CNVgc[i + 1]) / 2
    index = typeCNV != 0

    CNVRD = CNVRD[index]
    CNVchr = CNVchr[index]
    CNVstart = CNVstart[index]
    CNVend = CNVend[index]
    CNVlen = CNVlen[index]
    CNVmq = CNVmq[index]
    CNVgc = CNVgc[index]
    CNVtype = typeCNV[index]

    CNVdata = [*zip(CNVchr, CNVstart, CNVend, CNVlen, CNVtype)]
    final_CNV = pd.DataFrame(CNVdata, columns=['chr', 'start', 'end', 'size', 'type'])

    return final_CNV

'''
 use rd, gc and mq as features

'''


# get params

bam = sys.argv[1]
reference = sys.argv[2]
binSize = int(sys.argv[3])
reseg_len = int(sys.argv[4])
outfile = sys.argv[5]
alpha = 0.25

# get RD&MQ

chrList = get_chrlist(bam)
chrNum = len(chrList)
refList = [[] for i in range(22)]
refList,chr_num,chr_name = read_ref_file(reference, refList)
chrLen = np.full(22, 0)
for i in range(22):
    chrLen[i] = len(refList[i])

# Binning

print("Read bam file:", bam)
ReadCount = np.full(np.max(chrLen), 0)
Mapq = np.full(np.max(chrLen), 0)
ReadCount, bin_pos, Mapq = get_RC(bam, ReadCount, Mapq)
bin_start, bin_end, bin_len, bin_RD, bin_gc = ReadDepth(ReadCount, refList[chr_num-1], bin_pos)
print('mean_rd:',np.mean(bin_RD))

# CBS

with open(outfile + '.txt', 'w') as file:
    for c in range(len(bin_RD)):
        file.write(str(bin_RD[c]) + '\n')
subprocess.call('Rscript CBS_data.R ' + outfile, shell=True)
os.remove(outfile + '.txt')

# Segmentation

seg_start, seg_end, seg_len, reseg_gc = Read_seg_file(bin_start,bin_len,bin_end,bin_gc)
reseg_count, reseg_mq, reseg_start, reseg_end = resegment_RD(ReadCount, Mapq, seg_start, seg_end)
reseg_mq /= reseg_count
res_rd = prox_tv1d(alpha, reseg_count)
reseg_count = res_rd
reseg_chr = []
reseg_chr.extend(21 for j in range(len(reseg_count)))
data = [*zip(reseg_chr,reseg_start, reseg_end, reseg_count, reseg_gc, reseg_mq)]
data = pd.DataFrame(data, columns=['chr','start', 'end', 'rd', 'gc', 'mq'])

# PCC&call_CNVs

called_CNVs = PCC(data)
print('CNVs:' + '\n',called_CNVs)
with open(outfile + '.result.txt', 'w', ) as Outfile:
    called_CNVs.to_string(Outfile)
