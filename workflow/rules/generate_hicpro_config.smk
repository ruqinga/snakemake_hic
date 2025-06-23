# 允许参数未提供
# rule generate_hicpro_config:
#     input:
#         template="config/config-hicpro-template.j2"
#     output:
#         config_txt="config/config-hicpro.txt"
#     params:
#         config=config  # 将 Snakemake config 传给脚本
#     run:
#         from jinja2 import Template
#
#         with open(input.template) as f:
#             template = Template(f.read())
#
#         rendered = template.render(config=params.config)
#
#         with open(output.config_txt, "w") as out_f:
#             out_f.write(rendered)

# 严格模式，未提供参数会报错
rule generate_hicpro_config:
    input:
        template="config/config-hicpro-template.j2"
    output:
        config_txt="config/config-hicpro.txt"
    params:
        config=config
    run:
        from jinja2 import Template, StrictUndefined
        from jinja2.exceptions import UndefinedError

        try:
            with open(input.template) as f:
                template = Template(f.read(), undefined=StrictUndefined)

            rendered = template.render(config=params.config)

            with open(output.config_txt, "w") as out_f:
                out_f.write(rendered)

        except UndefinedError as e:
            print(f"[Error] Missing variable in config.yaml: {e}")
            raise e