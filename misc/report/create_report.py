import click
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from jinja2.exceptions import TemplateNotFound
import utils
import seaborn as sns


@click.command()
@click.argument("cf_results")
def create_report(cf_results):

    cf_results_df = pd.read_csv(cf_results, sep="\t")

    env = Environment(loader=FileSystemLoader("./"))
    template = env.get_template("template_report.html")

    template_vars = dict()

    variant_category(cf_results_df)
    dnm_inherited_unknown(cf_results_df)
    functional_consequences(cf_results_df)

    variants_per_gene(cf_results_df, template_vars)

    outputText = template.render(template_vars)
    with open("index.html", "w") as f:
        f.write(outputText)


def variant_category(cf_results_df):
    """
    Create an horizontal bar chart with the number of variants
    beeing SNV, indel or CNV.

    Args:
        cf_results_df (pd.DataFrame): clinical filtering results
    """

    nb_cnvs = (cf_results_df.cnv_length != ".").sum()
    nb_indels = (cf_results_df.indel_length != ".").sum()
    nb_snvs = cf_results_df.shape[0] - nb_cnvs - nb_indels

    assert nb_snvs == ((cf_results_df.indel_length == ".") & (cf_results_df.cnv_length == ".")).sum()

    data = {"Category": ["CNV", "Indel", "SNV"], "count": [nb_cnvs, nb_indels, nb_snvs]}
    df = pd.DataFrame(data)
    df.index = df.Category

    utils.horizontal_barchart(df, "Variant per category", "variant_category_distribution")


def dnm_inherited_unknown(cf_results_df):
    """
    Create an horizontal bar chart with the number of variants
    beeing inherited, DNM or unknown (singletons).

    Args:
        cf_results_df (pd.DataFrame): clinical filtering results
    """

    nb_dnms = cf_results_df.DNM.sum()
    nb_unknown = ((cf_results_df.mum == "0") & (cf_results_df.dad == "0")).sum()
    nb_inherited = cf_results_df.shape[0] - nb_dnms - nb_unknown

    # TODO : add an assertion on decipher_inheritance when fully ready

    data = {"Inheritance": ["DNM", "Inherited", "Unknown"], "count": [nb_dnms, nb_inherited, nb_unknown]}
    df = pd.DataFrame(data)
    df.index = df.Inheritance

    utils.horizontal_barchart(df, "Variant inheritance", "variant_inheritance_distribution")


def functional_consequences(cf_results_df):
    """
    Create an horizontal bar chart with the functional consequences for each variant category

    Args:
        cf_results_df (pd.DataFrame): clinical filtering results
    """

    cf_snv = cf_results_df.loc[(cf_results_df.cnv_length == ".") & (cf_results_df.indel_length == ".")]
    cf_indel = cf_results_df.loc[cf_results_df.indel_length != "."]
    cf_cnv = cf_results_df.loc[cf_results_df.cnv_length != "."]
    variant_category = ["SNV", "Indel", "CNV"]

    for idx, cf_df in enumerate([cf_snv, cf_indel, cf_cnv]):
        data = cf_df.consequence.value_counts()
        df = pd.DataFrame(data)
        df = df.sort_values(by="count", ascending=False)

        utils.horizontal_barchart(df, f"{variant_category[idx]}", f"{variant_category[idx]}_functional_consequence")


def variants_per_gene(cf_results_df, template_vars):
    """
    Prepare table that counts number of variants and number
    of DNM per gene

    Args:
        cf_results_df (pd.DataFrame): CF results
        template_vars (dict): Jinja2 template dictionnary

    """

    nb_variants_per_gene_df = pd.DataFrame(cf_results_df.symbol.value_counts())
    nb_variants_per_gene_df.columns = ["variants"]

    nb_dnm_per_gene_df = pd.DataFrame(cf_results_df.loc[cf_results_df.DNM == True].symbol.value_counts())
    nb_dnm_per_gene_df.columns = ["DNM"]

    df = nb_variants_per_gene_df.merge(nb_dnm_per_gene_df, on="symbol", how="outer")
    df.fillna(0, inplace=True)
    df["symbol"] = df.index
    df = df[["symbol", "variants", "DNM"]]
    df["DNM"] = df["DNM"].astype(int)

    html_table = utils.dfToHtmlTable(df)

    template_vars["table_variants_per_gene"] = html_table


if __name__ == "__main__":

    # Set the Seaborn style
    sns.set_style("whitegrid")

    create_report()
