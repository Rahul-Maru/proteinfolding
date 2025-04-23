import matplotlib.pyplot as plt
import logomaker as logomaker
import pandas as pd


def main():
    logomaker.list_font_names()

    repr_mat = pd.read_csv("hsm/outs/logo/repr.csv", index_col=0)
    repr_mat = repr_mat.fillna(0)

    # create Logo object
    crp_logo = logomaker.Logo(repr_mat,
                            shade_below=.5,
                            fade_below=.5,
                            font_name="Comic Sans MS")
    print(crp_logo.has_been_drawn)

    # style using Logo methods
    crp_logo.style_spines(visible=False)
    crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
    crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

    # style using Axes methods
    crp_logo.ax.set_ylabel("f", labelpad=-1)
    crp_logo.ax.xaxis.set_ticks_position('none')
    crp_logo.ax.xaxis.set_tick_params(pad=-1)

    # style and show figure
    crp_logo.fig.show()
    print(crp_logo.has_been_drawn)

    input("press any key to quit. ")


if __name__ == "__main__":
    main()
