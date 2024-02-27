import matplotlib.pyplot as plt
import seaborn as sns


def horizontal_barchart(df, title, filename):
    """
    Create an horizontal bar chart

    Args:
        df (pd.DataFrame): data
    """

    df = df.sort_values(by="count", ascending=False)

    # Create the horizontal bar plot
    plt.figure(figsize=(12, 6))
    ax = sns.barplot(x="count", y=df.index, data=df, color="skyblue")

    # Add labels and title with custom font settings
    plt.title(title, fontsize=16)

    # Add data values on the bars
    for index, value in enumerate(df["count"]):
        plt.text(value, index, str(value), va="center", fontsize=10, fontweight="bold")

    # Remove the top and right spines
    sns.despine()

    plt.tight_layout()
    plt.savefig(f"plots/{filename}.svg", format="svg")


def dfToHtmlTable(df):
    """
    Create the HTML code corresponding to a dataframe

    Args:
        df (pd.DataFrame): data

    Returns:
        str: The data frame as an HTML string
    """

    header = df.columns.values.tolist()

    html = '<table id="count_variants_genes" class="display" cellspacing="0" width="100%">'
    html += "<thead><tr>"
    for column in header:
        html += "<th>" + column + "</th>"
    html += "</tr></thead>"

    html += "<tfoot><tr>"
    for column in header:
        html += "<th>" + column + "</th>"
    html += "</tr></tfoot>"

    html += "<tbody>"
    for index, row in df.iterrows():
        rowList = row.tolist()
        html += "<tr>"
        for value in rowList:
            html += "<td>" + str(value) + "</td>"
        html += "</tr>"

    html += "</tbody>"
    html += "</table>"

    return html
