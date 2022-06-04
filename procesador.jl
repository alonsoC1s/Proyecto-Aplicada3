### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 0314c9ba-8e88-11ec-0b13-ff6d7516e4b5
begin
	using Pkg
	#Pkg.add(url="https://github.com/sdBrinkmann/HPFilter")
	Pkg.activate()
end

# ╔═╡ 84782bc2-92c5-4662-95f9-09c3312b6d4d
begin
	import GeoDataFrames as GDF
	using DataFramesMeta, CategoricalArrays
	using CSV, Statistics
	using LinearAlgebra
	using SparseArrays
	using SingularSpectrumAnalysis
	using Interpolations
	using NearestNeighbors
	using BenchmarkTools
	#using Plots
	using StatsPlots
	using Distances
	using ArchGDAL
	using Interpolations
end

# ╔═╡ 58cdca1b-0ae8-409c-ade4-272fc741b392
md"""
# Prototipo para ver si funciona el concepto de Aplicada 2
"""

# ╔═╡ b21f5c87-2700-4914-8ab4-8c1190e21148
md"## 1. Obteniendo, limpiando y acomodando datos"

# ╔═╡ 0644a198-509d-4423-a715-677cd65ee320
# ╠═╡ disabled = true
#=╠═╡
begin
	function segmenta_alcaldia(str_alcaldias)
		if occursin(r"Estado de México", str_alcaldias)
			return "EDOMX"
		elseif !occursin(r"/", str_alcaldias)
			return str_alcaldias
		else
			# Partiendo por "/" y regresando
			return split(str_alcaldias, "/")[rand(1:2)] |> rstrip |> lstrip
		end
	end

	
	cd("../DataLab/")
	mmetrobus = GDF.read("Metrobus_estaciones_utm14n.geojson")
	
	# Convirtiendo GEOM directo a lat, lng por facilidad
	# Convirtiendo a categorico lo que debería serlo
	mmetrobus = @chain mmetrobus begin
		@transform! @byrow :ALCALDIAS = segmenta_alcaldia(:ALCALDIAS)
		@transform! @byrow :SISTEMA = "Metrobús"
		@transform! begin
			:geom = ArchGDAL.setcoorddim!.(:geom, 2)
			:lng = ArchGDAL.getx.(:geom, 0)
			:lat = ArchGDAL.gety.(:geom, 0)
			$([:LINEA, :TIPO, :ALCALDIAS, :SISTEMA] .=> categorical
				.=> [:LINEA, :TIPO, :ALCALDIA, :SISTEMA]
			)
		end
		# Quitando accesibilidad para tener paridad de datos
		#@select! $(Not(:ACCESIB))
		@select! $[:NOMBRE, :ALCALDIA, :SISTEMA, :LINEA, :AÑO,:lng, :lat]
	end
	#CSV.write("../Proyecto-Aplicada2/data/metrobus_procesado.csv", mmetrobus; bom=true)
end
  ╠═╡ =#

# ╔═╡ 9a6ea0ba-8c89-4728-826c-2413ebf60b8a
#=╠═╡
 mmetrobus[1, :].geom
  ╠═╡ =#

# ╔═╡ 893be8c3-90e5-4498-bab5-2b500045f989
# ╠═╡ disabled = true
#=╠═╡
cd("../Proyecto-Aplicada2")
  ╠═╡ =#

# ╔═╡ 222349ef-d993-4af0-9975-f94acae44b88
md"Ahora tomo los datos del metrobus y lo junto todo"

# ╔═╡ cfc25d9d-35ae-48f8-83e5-bce6c415eea6
begin
	metro = CSV.read("data/metro_procesado.csv", DataFrame)
	metrobus = CSV.read("data/metrobus_procesado.csv", DataFrame)
	tren = CSV.read("data/trenligero_procesado.csv", DataFrame)
	cablebus = CSV.read("data/cablebus_procesado.csv", DataFrame)
end

# ╔═╡ 0ef2d15e-14ae-47cf-b73e-2565b621f670
# Juntando todo el transporte público
begin
	tpublico = vcat(metro, metrobus, cablebus)
	select!(tpublico, Not(:LINEA))
	# CSV.write("data/t_publico.csv", tpublico, bom=true)
end

# ╔═╡ 645ef6ea-4dd0-47f4-8624-d12654e383e4
begin
	const AÑOS = 1969:2021
	tpublico_dense = DataFrame(NOMBRE=String[],
		ALCALDIA=String[], SISTEMA=String[], AÑO=Int[],
		lng=Float64[], lat=Float64[]
	)
	
	for year in AÑOS
		# Obteniendo total de estaciones al año x
		tpu_yr = tpublico[tpublico.AÑO .<= year, :]
		tpu_yr.AÑO .= year
	
		append!(tpublico_dense, tpu_yr)
	end
end

# ╔═╡ 34ad2e05-1819-4a82-ae09-41b5c87743e6
tpublico_dense

# ╔═╡ 7b084347-1b81-4737-be99-88492903ee6b


# ╔═╡ d930feb1-d7dc-41df-8d53-478e2d1c5a71
md"Midiendo distancia de cada zona residencial a su estación de transporte público masivo más cercano por año y por alcaldía"

# ╔═╡ 0ff0079e-24eb-4edd-865c-8b9fd91215ea
uso_suelo = CSV.read("data/uso_de_suelo.csv", DataFrame)

# ╔═╡ 92866aa3-e980-4b60-b49b-cf2126edaecc
@subset(uso_suelo, @byrow startswith(:ALCALDIA, "Álvar"))

# ╔═╡ fb2f127d-c158-48c6-9762-583eed7bfbdf
md"#### Observación importante, no hay nada de datos para Álvaro Obregón"

# ╔═╡ bf354eac-7ed5-4713-828e-305de45f502a
md"
Continuo a medir la distancia promedio de las alcaldías que vienen en uso de suelo a su estación de transporte masivo más cercano

Estrategia, hago las operaciones sobre el DF de tpublico porque las zonas habitacionales son estáticas."

# ╔═╡ c28ac129-afd3-4071-b1ac-02d8b1c7565a
begin
	function nearest_neighbor(tree, points)
		_, dist = nn(tree, points)
		return dist
	end
	
	function distancia_promedio(tlats, tlngs)
		# Ajustando kdtree al transporte público existente
		kdtree = BallTree(hcat(tlats, tlngs)', Haversine(); leafsize=60)
	
		# Evaluando en las zonas habitacionales
		return @chain uso_suelo begin
			#_, dsts = nn(kdtree, adjoint(hcat(:lat, :lng)))
			@transform :distancia = nearest_neighbor(kdtree, hcat(:lat, :lng)')
			groupby(:ALCALDIA)
			@combine :dst_mean = mean(:distancia)
			permutedims(1) # Transponiendo df para volver observaciones columnas.
			@select $(Not(:ALCALDIA))
		end
	end
end

# ╔═╡ 5540ed04-e09d-4d56-8587-5b2afae5abc2
mean_dist_tp = @chain tpublico_dense begin
	# Primera transformación
	## Encontrando distancia promedio a transporte por alcaldía y año
	groupby(:AÑO)
	@combine $AsTable = distancia_promedio(:lat, :lng)
	stack()
	@transform $([:variable, :value] .=> [:ALCALDIA, :MEAN_DIST])
	@select $(Not([:variable, :value]))
	sort([:AÑO, :ALCALDIA])
end

# ╔═╡ 9375e07f-d0ee-412f-9822-635ce2e6772e
md"Ahora, contando el total de estaciones por alcaldía a cada año, y juntandolo con `datoz` con un inner join. Luego, poniendo población y alguna otra variable "

# ╔═╡ 6d256162-b01f-400e-a59b-033f5e1b96d4
# Ojo, aqui faltan las alcaldias que no tienen transporte e.g. Cuajimalpa
estaciones_p_año = @chain tpublico begin
	groupby([:AÑO, :ALCALDIA])
	@combine $(nrow => :EST_TOT)
end

# ╔═╡ 5641e8e1-8b86-4f94-b99c-c6cc63bd18be
# ╠═╡ show_logs = false
md" Para poder juntar con éxito con el resto de la info necesito la información para cada año intermedio entre construcciones nuevas.

Hago un dataframe vacio que contiene todos los años desde 1969 hasta 2021 con total de estaciones igual a cero. Después hago un innerjoin con año y alcaldía. Eso llena los espacios vacíos. Luego agrupo y reduzco con cumsum"

# ╔═╡ 1bbc9870-431f-4511-9292-62dd30ea1383
function densify_timeseries(df::DataFrame, col::Symbol)
	alcaldias = levels(df.ALCALDIA)

	len = length(AÑOS) * length(alcaldias)

	dense_alc = repeat(alcaldias; inner=length(AÑOS))
	dense_yrs = repeat(AÑOS; outer=length(alcaldias))

	denso = DataFrame(ALCALDIA=dense_alc, AÑO=dense_yrs, DUMMY = zeros(len))

	return @chain outerjoin(denso, df; on=[:ALCALDIA, :AÑO], makeunique=true) begin
		coalesce.(0)
		sort
		@transform! @byrow $col = $col + :DUMMY
		groupby(:ALCALDIA)
		@combine begin 
			:EST_TOTAL = cumsum($col)
			:AÑO
		end
	end
end

# ╔═╡ a38d84aa-5b3d-4fea-8c4e-aab2e5b4ea54
md"Midiendo distancia al centro de la ciudad"

# ╔═╡ 30046cd7-733e-4a38-8089-4092aa491976
begin
	zocalo = [19.432595, -99.133525]

	dist_zoc = @chain uso_suelo begin
		@transform @byrow :zdist = haversine(zocalo, [:lat, :lng])
		groupby(:ALCALDIA)
		@combine :ZOC_DIST = mean(:zdist)
		sort!(:ZOC_DIST)
	end
end

# ╔═╡ 2d7fe819-45b0-4094-a088-c6700fa8169c
md"""Una vez más tengo que "densificar" """

# ╔═╡ fd891304-043b-4441-a709-4113878d5514
function densify_timeseries(df::DataFrame)
	alcaldias = levels(df.ALCALDIA)

	len = length(AÑOS) * length(alcaldias)

	dense_alc = repeat(alcaldias; inner=length(AÑOS))
	dense_yrs = repeat(AÑOS; outer=length(alcaldias))
	values = repeat(df.ZOC_DIST; inner=length(AÑOS))

	return DataFrame(ALCALDIA=dense_alc, AÑO=dense_yrs, ZOC_DIST=values)
end

# ╔═╡ 55a35d06-de32-43f2-909e-1c6fc3ce6b3e
est_annum_dense = densify_timeseries(estaciones_p_año, :EST_TOT)

# ╔═╡ 6919245b-f9cb-42e6-87fe-1fab2067029d
# Verificando que el conteo coincide con el número de rows de tpublico
est_annum_dense[[53*i for i=1:14], :].EST_TOTAL |> sum

# ╔═╡ 8980ae50-3c0f-492e-a198-ecea34e1d348
dist_zoc_dense = densify_timeseries(dist_zoc)

# ╔═╡ 837d0c64-2b58-4f89-a305-2ddd88e5b3a5
md"""
Poniendo la población

El reto es la interpolación lineal y verificar que el join (outer?) no pierda años
"""

# ╔═╡ f332e25d-b0ca-46fa-92f3-e8de2c3c615d
begin
	cuajimalpizaciones = Dict(
		"Cuajimalpa de Morelos" => "Cuajimalpa",
		"La Magdalena Contreras" => "Magdalena Contreras",
	)
	
	function descuajimalpiza(alc)
		if alc in keys(cuajimalpizaciones)
			return cuajimalpizaciones[alc]
		else
			return alc
		end
	end

	# Leer población de 1970 y 1980
	pob_70_80 = CSV.read("data/pob70-80.csv", DataFrame)
	
	pob = @chain CSV.read("data/poblacion_alcaldias.csv", DataFrame) begin
		@subset @byrow :ALCALDIA != "CDMX"
		@transform begin
			@byrow :ALCALDIA = descuajimalpiza(:ALCALDIA)
		end
	end

	pob = vcat(pob, pob_70_80) |> sort
end

# ╔═╡ b58a28ab-3eaa-4a84-9485-a0f930082c7c
# Interpolación 
begin
	pob_groups = groupby(pob, :ALCALDIA)
	pob_llenito = DataFrame(AÑO=Int64[], ALCALDIA=String[], POBLACION=Float64[])
	
	for group_key in keys(pob_groups)
		pg = pob_groups[group_key]
		
		interp_linear = LinearInterpolation(
				pg.AÑO,
				pg.POBLACION;
				extrapolation_bc = Line()
		)
		años = 1969:1:2021
		
		dief = DataFrame(
			AÑO = años,
			ALCALDIA = group_key.ALCALDIA,
			POBLACION = interp_linear(años)
		)
		append!(pob_llenito, dief)
	end
end

# ╔═╡ c82c7bc6-c7cd-4376-b5ae-779a8d0ecb6e
pob_llenito

# ╔═╡ a868f876-453b-423c-ab61-72ca74f9099b
md"Juntando todo"

# ╔═╡ d827df20-73a0-4169-ab06-85ba302aff73
begin
	# Juntando todos los datos
	final = @chain begin
		leftjoin(pob_llenito, mean_dist_tp; on=[:AÑO, :ALCALDIA])
		leftjoin(est_annum_dense; on=[:AÑO, :ALCALDIA])
		leftjoin(dist_zoc_dense; on=[:AÑO, :ALCALDIA])
		# CSV.write("data/datosfinales.csv", _; bom=true)
	end
	coalesce(final.EST_TOTAL, 0)
end

# ╔═╡ a0d28b5e-7863-44ee-8abe-78fe1fae410f
final

# ╔═╡ 91a306e7-a4e4-438a-9e07-35e69ef5af92
# CSV.write("data/datosfinales.csv", final; bom=true)

# ╔═╡ 314da61a-6943-4f89-baf9-af140d795bd1
pwd()

# ╔═╡ Cell order:
# ╠═0314c9ba-8e88-11ec-0b13-ff6d7516e4b5
# ╠═84782bc2-92c5-4662-95f9-09c3312b6d4d
# ╠═58cdca1b-0ae8-409c-ade4-272fc741b392
# ╠═b21f5c87-2700-4914-8ab4-8c1190e21148
# ╠═0644a198-509d-4423-a715-677cd65ee320
# ╠═9a6ea0ba-8c89-4728-826c-2413ebf60b8a
# ╠═893be8c3-90e5-4498-bab5-2b500045f989
# ╟─222349ef-d993-4af0-9975-f94acae44b88
# ╠═cfc25d9d-35ae-48f8-83e5-bce6c415eea6
# ╠═0ef2d15e-14ae-47cf-b73e-2565b621f670
# ╠═645ef6ea-4dd0-47f4-8624-d12654e383e4
# ╠═34ad2e05-1819-4a82-ae09-41b5c87743e6
# ╠═7b084347-1b81-4737-be99-88492903ee6b
# ╟─d930feb1-d7dc-41df-8d53-478e2d1c5a71
# ╠═0ff0079e-24eb-4edd-865c-8b9fd91215ea
# ╠═92866aa3-e980-4b60-b49b-cf2126edaecc
# ╟─fb2f127d-c158-48c6-9762-583eed7bfbdf
# ╟─bf354eac-7ed5-4713-828e-305de45f502a
# ╠═c28ac129-afd3-4071-b1ac-02d8b1c7565a
# ╠═5540ed04-e09d-4d56-8587-5b2afae5abc2
# ╟─9375e07f-d0ee-412f-9822-635ce2e6772e
# ╠═6d256162-b01f-400e-a59b-033f5e1b96d4
# ╟─5641e8e1-8b86-4f94-b99c-c6cc63bd18be
# ╠═1bbc9870-431f-4511-9292-62dd30ea1383
# ╠═55a35d06-de32-43f2-909e-1c6fc3ce6b3e
# ╠═6919245b-f9cb-42e6-87fe-1fab2067029d
# ╟─a38d84aa-5b3d-4fea-8c4e-aab2e5b4ea54
# ╠═30046cd7-733e-4a38-8089-4092aa491976
# ╟─2d7fe819-45b0-4094-a088-c6700fa8169c
# ╠═fd891304-043b-4441-a709-4113878d5514
# ╠═8980ae50-3c0f-492e-a198-ecea34e1d348
# ╟─837d0c64-2b58-4f89-a305-2ddd88e5b3a5
# ╠═f332e25d-b0ca-46fa-92f3-e8de2c3c615d
# ╠═b58a28ab-3eaa-4a84-9485-a0f930082c7c
# ╠═c82c7bc6-c7cd-4376-b5ae-779a8d0ecb6e
# ╟─a868f876-453b-423c-ab61-72ca74f9099b
# ╠═d827df20-73a0-4169-ab06-85ba302aff73
# ╠═a0d28b5e-7863-44ee-8abe-78fe1fae410f
# ╠═91a306e7-a4e4-438a-9e07-35e69ef5af92
# ╠═314da61a-6943-4f89-baf9-af140d795bd1
