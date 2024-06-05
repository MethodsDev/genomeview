from jinja2 import Environment, PackageLoader, select_autoescape



def get_jinja_environment():
    env = Environment(
        loader=PackageLoader('genomeview', 'templates'),
        autoescape=select_autoescape(['html', 'xml'])
    )
    return env



def render_tab_titles(tabs, page_title):
    env = get_jinja_environment()
    template = env.get_template('tab_titles.html')
    return template.render(tabs=tabs, page_title=page_title)

def render_bam_buttons_views(tab):
    env = get_jinja_environment()
    template = env.get_template('bam_buttons_views.html')
    return template.render(tab=tab)

def render_bam_views(bams):
    env = get_jinja_environment()
    template = env.get_template('bam_views.html')
    return template.render(bams=bams)


def plot_sorted_support_as_tabs(config, bams_dict, page_title, **kwargs):
    tabs = []
    for feature_name, virtual_bam_dict in bams_dict.items():
        interval = config.id_to_coordinates[feature_name]
        feature_name = config.gene_id_to_gene_name[config.transcript_to_gene[feature_name]] + "_" + feature_name

        bams = []
        for classification, virtual_bam in virtual_bam_dict.items():
            if virtual_bam is None:
                continue
            unique_id = f"{feature_name}_{classification}"
            coverage_svg = config.plot_interval(
                bams_dict={classification: virtual_bam},
                interval=interval,
                with_reads=False,
                with_coverage=True,
                add_track_label=False,
                fill_coverage=True,
                **kwargs
            )._repr_svg_()

            reads_svg = config.plot_interval(
                bams_dict={classification: virtual_bam},
                interval=interval,
                with_reads=True,
                with_coverage=False,
                with_axis=False,
                with_bed=False,
                add_track_label=False,
                add_reads_label=False,
                vertical_layout_reads=True,
                **kwargs
            )._repr_svg_()

            bams.append({
                'unique_id': unique_id,
                'name': f"{feature_name}_{classification}",
                'statis_svg': coverage_svg,
                'resizable_svg': reads_svg
            })

        tabs.append({
            'feature_name': feature_name,
            'bams': bams
        })

    return render_tab_titles(tabs, page_title)
