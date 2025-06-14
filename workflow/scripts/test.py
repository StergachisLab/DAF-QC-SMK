regions = snakemake.params.regions



with open(snakemake.output[0], 'w') as out_f:
	out_f.write(f"Regions: {regions}\n")
	out_f.write(f"Type of regions: {type(regions)}\n")
	regions.split(',')
	for region in regions.split(','):
		out_f.write(f"Region: {region}\n")
		region_parts = region.replace(':', '-').split('-')
		out_f.write(f"Region parts: {region_parts}\n")
	

