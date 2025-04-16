import matplotlib.pyplot as plt
import logomaker as logomaker


def main():
    logomaker.list_font_names()

    # load crp energy matrix
    crp_df = -logomaker.get_example_matrix('crp_energy_matrix',
                                            print_description=False)

    print((crp_df))
    # create Logo object
    crp_logo = logomaker.Logo(crp_df,
                            shade_below=.5,
                            fade_below=.5,
                            font_name="Comic Sans MS")
    print(crp_logo.has_been_drawn)

    # style using Logo methods
    # crp_logo.style_spines(visible=False)
    # crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
    # crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

    # # style using Axes methods
    # crp_logo.ax.set_ylabel(r"$-\Delta \Delta G$ (kcal/mol)", labelpad=-1)
    # crp_logo.ax.xaxis.set_ticks_position('none')
    # crp_logo.ax.xaxis.set_tick_params(pad=-1)

    # style and show figure
    crp_logo.fig.show()
    print(crp_logo.has_been_drawn)

    input("press any key to quit")


if __name__ == "__main__":
    main()
