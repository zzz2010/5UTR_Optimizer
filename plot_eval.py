from matplotlib import pyplot as plt
import pandas as pd
import sys,os
import seaborn as sns

if __name__ == '__main__':
    excel_fn=sys.argv[1]
    TE_result_df=pd.read_excel(excel_fn,sheet_name="TE_result")

    RNA_result_df = pd.read_excel(excel_fn, sheet_name="RNA_result")
    fig,axes=plt.subplots(nrows=2,ncols=2,figsize=(18,20))
    axes=axes.flatten()
    palette="Set1"
    ax=sns.boxplot(x="model", y="sp_cor",
                hue="model", palette=palette,
                data=TE_result_df,ax=axes[0])
    ax.set_title("TE Prediction Spearman Correlation (10-fold CVs)")
    ax.set_ylabel("TE Prediction Spearman Correlation")

    ax=sns.boxplot(x="model", y="rsq",
                hue="model", palette=palette,
                data=TE_result_df,ax=axes[1])
    ax.set_title("TE Prediction R^2 (10-fold CVs)")
    ax.set_ylabel("TE Prediction R^2")


    ax=sns.boxplot(x="model", y="sp_cor",
                hue="model", palette=palette,
                data=RNA_result_df,ax=axes[2])
    ax.set_title("RNA RPKM Prediction Spearman Correlation (10-fold CVs)")
    ax.set_ylabel("RNA RPKMn Spearman Correlation")

    ax=sns.boxplot(x="model", y="rsq",
                hue="model", palette=palette,
                data=RNA_result_df,ax=axes[3])
    ax.set_title("RNA RPKM Prediction R^2 (10-fold CVs)")
    ax.set_ylabel("RNA RPKM  Prediction R^2")

    plt.suptitle(os.path.basename(excel_fn))
    plt.show()