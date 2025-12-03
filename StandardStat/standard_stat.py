import sys
import pandas as pd


def main():
    concordance = pd.read_csv(sys.argv[1], sep='\t')
    sample = sys.argv[2]
    bed = pd.read_csv(sys.argv[3], sep='\t', header=None)
    tp = concordance['TP'].sum()
    fp = concordance['FP'].sum()
    fn = concordance['FN'].sum()
    n = (bed[2] - bed[1]).sum() - tp - fn
    tn = n - fp
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    tp_snp = concordance.loc[concordance['type'] == 'SNP', 'TP'].values[0]
    tp_indel = tp - tp_snp
    fp_snp = concordance.loc[concordance['type'] == 'SNP', 'FP'].values[0]
    fp_indel = fp - fp_snp
    fn_snp = concordance.loc[concordance['type'] == 'SNP', 'FN'].values[0]
    fn_indel = fn - fn_snp
    tn_snp = tn - fp_snp
    tn_indel = tn - fp_indel
    sensitivity_snp = tp_snp / (tp_snp + fn_snp)
    sensitivity_indel = tp_indel / (tp_indel + fn_indel)
    specificity_snp = tn_snp / (tn_snp + fp_snp)
    specificity_indel = tn_indel / (tn_indel + fp_indel)
    accuracy_snp = (tp_snp + tn_snp) / (tp_snp + tn_snp + fp_snp + fn_snp)
    accuracy_indel = (tp_indel + tn_indel) / (tp_indel + tn_indel + fp_indel + fn_indel)
    dic = {'sample_id': sample, 'sensitivity': sensitivity, 'specificity': specificity, 'accuracy': accuracy, 'sensitivity_snp': sensitivity_snp, 'sensitivity_indel': sensitivity_indel, 'specificity_snp': specificity_snp, 'specificity_indel': specificity_indel, 'accuracy_snp': accuracy_snp, 'accuracy_indel': accuracy_indel, 'tp_snp': tp_snp, 'tp_indel': tp_indel, 'fp_snp': fp_snp, 'fp_indel': fp_indel, 'fn_snp': fn_snp, 'fn_indel': fn_indel}
    df = pd.DataFrame([dic])
    df.to_json(sys.argv[4], orient='records', force_ascii=False, lines=True)


if __name__ == '__main__':
    main()

