### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 4e2986be-e8df-11ec-1371-296feeb95fca
begin
	using Pkg
	Pkg.activate()
end

# ╔═╡ 6ce9af88-4ed5-4773-8dff-2fd3bb6b83c5
begin
	using Plots, ArchGDAL, DataFrames, DataFramesMeta, CSV, StatsPlots
	import GeoDataFrames as GDF
	theme(:ggplot2)
end

# ╔═╡ d6332502-567b-434e-aa83-2bede0b4ba16
pwd()

# ╔═╡ f42840bf-3dfb-4927-9183-03ca7e6b0abb
begin
	alc = GDF.read("data/alcaldias_cdmx/alcaldías_cdmx/alcaldias_cdmx.shp")
	
	renames = Dict(
	    "MILPA ALTA" => "Milpa Alta",
	    "BENITO JUÁREZ" => "Benito Juárez",
	    "GUSTAVO A. MADERO" => "Gustavo A. Madero",
	    "COYOACÁN" => "Coyoacán",
	    "MIGUEL HIDALGO" => "Miguel Hidalgo",
	    "LA MAGDALENA CONTRERAS" => "Magdalena Contreras",
	    "TLÁHUAC" => "Tláhuac",
	    "AZCAPOTZALCO" => "Azcapotzalco",
	    "IZTACALCO" => "Iztacalco", 
	    "ÁLVARO OBREGÓN" => "Álvaro Obregón",
	    "XOCHIMILCO" => "Xochimilco",
	    "VENUSTIANO CARRANZA" => "Venustiano Carranza",
	    "TLALPAN" => "Tlalpan",
	    "CUAJIMALPA DE MORELOS" => "Cuajimalpa",
	    "CUAUHTÉMOC" => "Cuauhtémoc", 
	    "IZTAPALAPA" => "Iztapalapa"
	)
	
	@transform! alc @byrow :ALCALDIA = renames[:nomgeo]

	mb = GDF.read("data/mb_shp/Metrobus_lineas_rutas_LT.shp")
	metro = GDF.read("data/mb_shp/Metrobus_lineas_rutas_LT.shp")
	tlig = GDF.read("data/ste_tren_ligero_shp/STE_TrenLigero_linea_utm14n.shp")
	cbus = GDF.read("data/ste_cablebus_shp/STE_Cablebus_lineas_utm14n.shp")
end

# ╔═╡ d143f821-5419-476a-967e-56b3abf14482
begin
	plot(alc.geom, fill="lightgray")
end

# ╔═╡ b3ba7669-0f92-4c01-87d8-91f4fc22290f
info = CSV.read("data/datosfinales.csv", DataFrame)

# ╔═╡ 50e88957-6548-48ff-88b8-2c0e756529eb
md"## Alcaldía por población"

# ╔═╡ 3db00b7a-5673-4e03-b5c9-a1fe3a16f581
begin
	C(g::ColorGradient) = RGB[g[z] for z=LinRange(0,1,30)]
	g = :inferno

	cgrad([colorant"#132B43", colorant"#56B1F7"])
	theme(:ggplot2; c = cgrad([colorant"#132B43", colorant"#56B1F7"], scale=:log))
end

# ╔═╡ 045c7228-b956-46fd-8f25-b19183678da2
@chain info begin
	@subset @byrow :AÑO == 2021
	#groupby(:ALCALDIA)
	rightjoin(alc, on=:ALCALDIA)
	@df _ plot(:geom, fill_z=:POBLACION')
end

# ╔═╡ 3cb1f784-eccb-49c1-be00-ccdfb5569604
md"Por año"

# ╔═╡ 215b149c-4b70-4750-b873-ec53779be45c
begin
	anim = @animate for year in levels(aug.AÑO)
		@chain info begin
			@subset @byrow :AÑO == year
			rightjoin(alc, on=:ALCALDIA)
			@df _ plot(:geom, fill_z=:POBLACION',
				title = "Alcaldías por población al año $(year)"
			)
		end
	end
	
	gif(anim, "presentacion_files/figure-revealjs/mapapob.gif"; fps=2)
end

# ╔═╡ cd095986-7a9b-4e38-baec-a16d59d88b64
begin
	anim2 = @animate for year in levels(aug.AÑO)
		@chain aug begin
			@subset @byrow :AÑO == year
			@subset @byrow :EST_TOTAL != 0
			@df _ plot(:geom, fill_z=:EST_TOTAL',
				title = "Estaciones por alcaldía al año $(year)",
				xlims=(-99.4, -98.9),
				ylims=(19, 19.6),
			)
		end
	end
	
	gif(anim2, "presentacion_files/figure-revealjs/mapaest.gif"; fps=2)
end

# ╔═╡ de9f53d8-4b4a-4756-84f9-3708bb35dfd2
begin
	anim3 = @animate for year in levels(aug.AÑO)
		@chain aug begin
			@subset @byrow :AÑO == year
			#groupby(:ALCALDIA)
			#rightjoin(alc, on=:ALCALDIA)
			#@subset @byrow :EST_TOTAL != 0
			@df dropmissing(_) plot(:geom, fill_z=:MEAN_DIST',
				title = "Distancia promedio al transporte al $(year)",
				xlims=(-99.4, -98.9),
				ylims=(19, 19.6),
				rev=true
			)
		end
	end
	
	gif(anim3, "presentacion_files/figure-revealjs/mapadist.gif"; fps=2)
end

# ╔═╡ 01e3cf4e-59b9-4ed5-b6fe-164675158ec6
md"### conectividad"

# ╔═╡ bf7b4d8b-5b1b-46fa-97ae-ec2d0ffa9c93
con_by_yr = CSV.read("data/con_by_year.csv", DataFrame)

# ╔═╡ 4204ea92-70cf-43c3-8193-b25f745eed91
aug_score = leftjoin(con_by_yr, alc; on=:ALCALDIA)

# ╔═╡ d9a94d0d-7d48-4f07-bcba-ea290e062eb9
begin
	anim4 = @animate for year in levels(aug_score.AÑO)
		@chain aug_score begin
			@subset @byrow :AÑO == year
			#groupby(:ALCALDIA)
			#rightjoin(alc, on=:ALCALDIA)
			#@subset @byrow :EST_TOTAL != 0
			@df _ plot(:geom, fill_z=:score',
				title = "Conectividad al $(year)",
			)
		end
	end
	
	gif(anim4, "presentacion_files/figure-revealjs/mapaconn.gif"; fps=2)
end

# ╔═╡ cd52b1ef-cfcb-407d-92cf-708f6cbedb18
begin
	theme(:ggplot2)
	
	anim5 = @animate for year in levels(aug_score.AÑO)
		@chain aug_score begin
				@subset @byrow :AÑO <= year
				@df _ plot(:AÑO, :score,
					group = :ALCALDIA,
					legend = :outertopright,
					xlims = (1969, 2021),
					title = "Conectividad al $(year)",
				)
		end
	end

	gif(anim5, "presentacion_files/figure-revealjs/connline.gif"; fps=2)
end

# ╔═╡ 79985177-3bcb-467b-a8db-d6e464a7a901


# ╔═╡ Cell order:
# ╠═4e2986be-e8df-11ec-1371-296feeb95fca
# ╠═6ce9af88-4ed5-4773-8dff-2fd3bb6b83c5
# ╠═d6332502-567b-434e-aa83-2bede0b4ba16
# ╠═f42840bf-3dfb-4927-9183-03ca7e6b0abb
# ╠═d143f821-5419-476a-967e-56b3abf14482
# ╠═b3ba7669-0f92-4c01-87d8-91f4fc22290f
# ╠═50e88957-6548-48ff-88b8-2c0e756529eb
# ╠═3db00b7a-5673-4e03-b5c9-a1fe3a16f581
# ╠═045c7228-b956-46fd-8f25-b19183678da2
# ╠═3cb1f784-eccb-49c1-be00-ccdfb5569604
# ╠═215b149c-4b70-4750-b873-ec53779be45c
# ╠═cd095986-7a9b-4e38-baec-a16d59d88b64
# ╠═de9f53d8-4b4a-4756-84f9-3708bb35dfd2
# ╠═01e3cf4e-59b9-4ed5-b6fe-164675158ec6
# ╠═bf7b4d8b-5b1b-46fa-97ae-ec2d0ffa9c93
# ╠═4204ea92-70cf-43c3-8193-b25f745eed91
# ╠═d9a94d0d-7d48-4f07-bcba-ea290e062eb9
# ╠═cd52b1ef-cfcb-407d-92cf-708f6cbedb18
# ╠═79985177-3bcb-467b-a8db-d6e464a7a901
