import sys
import cyvcf2
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from statannot import add_stat_annotation


def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description="This script was " + 
            "designed to plot the distribution of a given phenotype " +
            "for all genotypes in a VCF given a specific chromosome " + 
            "and position.\n\n")

    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument(
            "--gzvcf", 
            type=str, 
            required=True, 
            help="The name of the input gzipped VCF file. Must be tbi indexed.",
            action="store")

    requiredNamed.add_argument(
            "--phenotype-file",
            type=str,
            required=True,
            help="A CSV of relevant phenotypes with taxa in the first column." + 
            " Note: headers are **required**.",
            action="store")

    requiredNamed.add_argument(
            "--phenotype",
            type=str,
            required=True,
            help="The phenotype column to use in the phenotype file",
            action="store")

    requiredNamed.add_argument(
            "--position",
            type=str,
            required=True,
            help="Position of interest in Chr01_Position format (or what matches the VCF).",
            action="store")

    requiredNamed.add_argument(
            "--output",
            type=str,
            required=True,
            help="Name of the output figure",
            action="store")

    parser.add_argument(
            "--test",
            type=str,
            required=False,
            default="t-test_ind",
            help="Test type: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, " + 
            "Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal.",
            action="store")

    parser.add_argument(
            "--point-size",
            type=int,
            required=False,
            default=2,
            help="Size for points in swarmplot. If size=0, then points are removed.",
            action="store")

    parser.add_argument(
            "--groups",
            type=str,
            required=False,
            default=None,
            help="A CSV containing taxa,groupid.",
            action="store")

    parser.add_argument(
            "--highlights",
            type=str,
            required=False,
            default="all",
            help="A comma-separated list of groups to highlight. Ex. durra,milo",
            action="store")

    parser.add_argument(
            "--geneid",
            type=str,
            required=False,
            default="Sbicolor v3",
            help="Name of gene to add to title.",
            action="store")

    parser.add_argument(
            "--sample-of-interest",
            type=str,
            required=False,
            default=None,
            help="File containing samples of interest to print genotype information. Header (CommonName,CUSO,PI).",
            action="store")

    return parser.parse_args()


def get_df_from_VCF(filename, chrom, pos):
    data_matrix = []
    vcf_reader = cyvcf2.Reader(filename, mode=u'rb', threads=4)
    vcf_reader.set_index(filename + ".tbi")
    chrom_header = [i for i in vcf_reader.raw_header.split("\n") if i.startswith("#CHROM")][0].split()
    data_matrix.append(chrom_header)
    for record in vcf_reader(f"{chrom}:{pos}-{int(pos)+1}"):
        data_matrix.append(str(record).rstrip().split())
    df = pd.DataFrame(data_matrix[1:], columns=data_matrix[0])
    df = df.applymap(lambda x: x.split(":")[0].replace("|","/"))
    return df


def get_plot_df(df):
    locus = df[(df["#CHROM"]==chrom) & (df["POS"]==pos)]
    plot_df = locus.transpose().merge(phenotypes[phenotype_abbrv], left_index=True , right_index=True)
    plot_df.columns = ["Genotype", phenotype_abbrv]
    return plot_df


if __name__ == "__main__":
    args = parse_arguments()
    chrom, pos = args.position.split("_")
    point_size = args.point_size

    # Phenotype
    phenotype_abbrv = args.phenotype
    phenotypes = pd.read_csv(args.phenotype_file, index_col=0)
    phenotypes.head()
    
    # VCF
    df = get_df_from_VCF(args.gzvcf, chrom, pos)
    plot_df = get_plot_df(df) 
    # Drop missing
    plot_df = plot_df.loc[plot_df["Genotype"] != "./.",]
    # Drop hets
    # plot_df = plot_df.loc[plot_df["Genotype"] != "0/1",]
   
    order = sorted(plot_df.Genotype.unique())
    pairwise_tests = list(combinations(order, 2))
    if args.groups:
        groups = pd.read_csv(args.groups, header=None, index_col=0)
        groups.columns = ["Subpopulation"]
        if args.highlights != "all":
            highlight_groups = args.highlights.split(",")
            groups = groups[groups.Subpopulation.isin(highlight_groups)]
            plot_df = plot_df.merge(groups, left_index=True, right_index=True, how="outer")
            plot_df.Subpopulation = plot_df.Subpopulation.fillna("other")
            sns.swarmplot(x="Genotype", y=phenotype_abbrv, data=plot_df, 
                order=order, size=point_size, hue="Subpopulation", hue_order=["other"] + highlight_groups)
        else:
            plot_df = plot_df.merge(groups, left_index=True, right_index=True)
            sns.swarmplot(x="Genotype", y=phenotype_abbrv, data=plot_df,
                order=order, size=point_size, hue="Subpopulation")
    else:
        sns.swarmplot(x="Genotype", y=phenotype_abbrv, data=plot_df, order=order, size=point_size)
    ax = sns.boxplot(x="Genotype", y=phenotype_abbrv, data=plot_df, order=order, color='white')
    add_stat_annotation(ax, data=plot_df, x="Genotype", y=phenotype_abbrv, order=order,
            box_pairs=pairwise_tests, test=args.test, text_format='star', loc='outside', verbose=2)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.xlabel(args.geneid + " " + args.position)
    plt.tight_layout()
    plt.savefig(args.output)
    plt.close()

    plot_df.to_csv(args.output.split('.')[0] + ".plot_df.csv")    
    if args.sample_of_interest:
        samples_of_interest = pd.read_csv(args.sample_of_interest, index_col=1)
        for sample in list(samples_of_interest.index):
            print(samples_of_interest.loc[sample, "CommonName"] + " " + plot_df.loc[sample, "Genotype"])
