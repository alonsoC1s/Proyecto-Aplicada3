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
	using CSV, Glob, Dates, Statistics, StatsPlots
	using LinearAlgebra
	using SparseArrays
	using SingularSpectrumAnalysis
	using Interpolations
	using NearestNeighbors
end

# ╔═╡ a27b24b4-e9a5-419e-bcea-25b8cfbffc66
begin
	using BenchmarkTools
	data = rand(3, 10^4)
	k = 3
	point = rand(3)
	
	btr = BruteTree(data)
	ktr = KDTree(data)
	batr = BallTree(data)
	
	@benchmark idxs, dists = nn($btr, $point)
end

# ╔═╡ 58cdca1b-0ae8-409c-ade4-272fc741b392
md"""
# Prototipo para ver si funciona el concepto de Aplicada 2

**Objetivo**:
- Ver si hay correlaciones entre las variables: Número de estaciones de transporte público & utilización vehicular (o un proxy de eso).

Estoy usando solo metro y metrobus pero podemos tener acceso a todo [Afluencia Preliminar en Transporte Publico](https://datos.cdmx.gob.mx/dataset/afluencia-preliminar-en-transporte-publico/resource/c2a29b87-0809-4448-90a2-158d8fc5d180) & también tenemos datos de [Movilidad](https://datos.cdmx.gob.mx/dataset/movilidad-historico-covid-19) (específicamente para el contexto COVID-19 y cambios porcentuales).

### Pregunta
Se puede establecer una relación estilo: más transporte público $\implies$ menos viajes en automóvil?

- A nivel alcaldía
"""

# ╔═╡ b21f5c87-2700-4914-8ab4-8c1190e21148
md"## 1. Obteniendo, limpiando y acomodando datos"

# ╔═╡ 222349ef-d993-4af0-9975-f94acae44b88
md"Ahora tomo los datos del metrobus y lo junto todo"

# ╔═╡ 893be8c3-90e5-4498-bab5-2b500045f989
cd("../Proyecto-Aplicada2")

# ╔═╡ cfc25d9d-35ae-48f8-83e5-bce6c415eea6
begin
	metro = CSV.read("data/metro_procesado.csv", DataFrame)
	metrobus = CSV.read("data/metrobus_procesado.csv", DataFrame)
	tren = CSV.read("data/trenligero_procesado.csv", DataFrame)
	cablebus = CSV.read("data/cablebus_procesado.csv", DataFrame)
end

# ╔═╡ 0ef2d15e-14ae-47cf-b73e-2565b621f670
# Juntando todo electrico
begin
	tpublico = vcat(metro, tren, metrobus, cablebus)
	# CSV.write("data/t_publico.csv", tpublico, bom=true)
end

# ╔═╡ 32d95ebc-ae09-4cec-ab64-2e495032a008
groupby(tpublico, [:ALCALDIA, :AÑO])

# ╔═╡ d930feb1-d7dc-41df-8d53-478e2d1c5a71
md"Midiendo distancia de cada zona residencial a su estación de transporte público masivo más cercano por año y por alcaldía"

# ╔═╡ 0ff0079e-24eb-4edd-865c-8b9fd91215ea
uso_suelo = CSV.read("data/uso_de_suelo.csv", DataFrame)

# ╔═╡ 5110c22b-cebe-460f-9a0a-b65f2643c866
begin
	tg = groupby(tpublico, [:AÑO, :ALCALDIA])
	ug = groupby(uso_suelo, :ALCALDIA)

	ktg = keys(tg)
	ktg1 = ktg[1]
	tg[ktg1]
end

# ╔═╡ bf354eac-7ed5-4713-828e-305de45f502a
md"Estrategia, hago las operaciones sobre el DF de tpublico porque las zonas habitacionales son estáticas. El reto será encontrar una manera ergonómica de agrupar por alcaldía e indexar por esa misma alcaldía para hacer el cálculo sobre dos dataframes disjuntos."

# ╔═╡ c28ac129-afd3-4071-b1ac-02d8b1c7565a
function distancia_promedio(alc, tlats, tlngs)
	# Ajustando kdtree a los metros existentes
	kdtree = KDTree(hcat(tlats, tlngs), leafsize=10)

	# Evaluando en las zonas habitacionales
	_, dsts = nn(kdtree, hcat(uso_suelo.lat, uso_suelo.lng))
	#=
	@chain uso_suelo begin
		_, dsts = nn(kdtree, hcat(:lat, :lng))
		@transform :distancia = dsts
		#=
		groupby(:ALCALDIA)
		@combine :dst_mean = mean(:distancia)
		=#
	end
	=#
end

# ╔═╡ 04a38ace-8c57-401a-8759-43d9d77175d0


# ╔═╡ d35505d9-2fcd-43a7-9b06-b5d26be1e349
@benchmark idxs, dists = nn($ktr, $point)

# ╔═╡ 491f2be2-799f-42a5-8de7-d7fed574d3ce
@benchmark idxs, dists = nn($batr, $point)

# ╔═╡ 5540ed04-e09d-4d56-8587-5b2afae5abc2
@transform tg[ktg1] begin
	:tst = distancia_promedio(:ALCALDIA, :lat, :lng)
end

# ╔═╡ 5641e8e1-8b86-4f94-b99c-c6cc63bd18be
md"Paso de la muerte. Hago un dataframe vacio que contiene todos los años desde 1980 hasta 2022 con total de estaciones verdes igual a cero. Después hago un innerjoin con año y alcaldía. Eso llena los espacios vacíos. Luego agrupo y reduzco con cumsum"

# ╔═╡ 1bbc9870-431f-4511-9292-62dd30ea1383
function densify_timeseries(df::DataFrame, col::Symbol)
	años = 1980:2022
	alcaldias = levels(df.ALCALDIA)

	len = length(años) * length(alcaldias)

	dense_alc = repeat(alcaldias; inner=length(años))
	dense_yrs = repeat(años; outer=length(alcaldias))

	denso = DataFrame(ALCALDIA=dense_alc, AÑO=dense_yrs, DUMMY = zeros(len))

	return @chain outerjoin(denso, df; on=[:ALCALDIA, :AÑO], makeunique=true) begin
		coalesce.(0)
		sort
		@transform! col = :col + :DUMMY
		groupby(:ALCALDIA)
		@combine begin 
			:EST_TOTAL = cumsum(:EST_V)
			:AÑO
		end
	end
end

# ╔═╡ 4a071190-96ba-48c3-a262-5487f4c6155d
densify_timeseries(tpublico)

# ╔═╡ dc4208eb-e975-4c46-9d06-419e20bd1505
begin
	años = 1980:2022
	#alcaldias = levels(telectrico.ALCALDIA)

	#len = length(años) * length(alcaldias)

	#dense_alc = repeat(alcaldias; inner=length(años))
	#dense_yrs = repeat(años; outer=length(alcaldias))

	#denso = DataFrame(ALCALDIA=dense_alc, AÑO=dense_yrs, EST_V = zeros(len))
	denso = DataFrame(AÑO=años, EST_V = zeros(length(años)))
end

# ╔═╡ af5e5030-fbc4-4e7f-8b90-8fe214892290
#=
telectrico_denso = @chain outerjoin(denso, telectrico_chico; on=[:ALCALDIA, :AÑO], makeunique=true) begin
	coalesce.(0)
	sort
	@transform! :EST_V = :EST_V + :EST_V_1
	groupby(:ALCALDIA)
	@combine begin 
		:EST_V = cumsum(:EST_V)
		:AÑO
	end
end
=#

# ╔═╡ 0cb4cac9-d546-4f1f-ac46-9db51b990baa
telectricocdmx_denso = @chain outerjoin(denso, telectrico_chico; on=:AÑO, makeunique=true) begin
	coalesce.(0)
	sort
	@transform! :EST_V = :EST_V + :EST_V_1
	# groupby(:ALCALDIA)
	@combine begin 
		:EST_V = cumsum(:EST_V)
		:AÑO
	end
end

# ╔═╡ 993cb8bf-aba6-4b3e-888c-d3c639967a52
md"## Población"

# ╔═╡ 438104ae-88bd-45b9-b87d-9a16938d35c1
poblacion = @chain CSV.read("poblacion_alcaldias.csv", DataFrame) begin
	# @subset @byrow :ALCALDIA == "CDMX"
end

# ╔═╡ fa1bbb87-ed47-4ac5-b9b5-7cd511bd8349
begin
	#@df poblacion plot(:AÑO, :POBLACION, group=:ALCALDIA)
	#@df poblacion plot(LinearInterpolation(:AÑO, :POBLACION))
end

# ╔═╡ 71209b46-ecf2-4721-9e04-f64e76411b64
# ╠═╡ disabled = true
#=╠═╡
begin
	itp = interpolate((poblacion.AÑO,), poblacion.POBLACION, Gridded(Linear()))
	
	plot!(1990:1:2020, itp(1990:1:2020))
end
  ╠═╡ =#

# ╔═╡ 6ddf350b-b847-4867-afd5-d7bdab5fb62d
# ╠═╡ disabled = true
#=╠═╡
poblacion_denso = DataFrame(
	AÑO=1990:1:2020, ALCALDIA=fill("CDMX", 31),
	POBLACION=Float64.(itp(1990:1:2020))
)
  ╠═╡ =#

# ╔═╡ 3f725649-ae2c-41c2-899e-b4bcb9ee273b
# CSV.write("parque_1980-2020.csv", vehic; bom=true)

# ╔═╡ 15060361-9c2f-4557-b45c-1e893381d398
# ╠═╡ disabled = true
#=╠═╡
begin
	MAIN = innerjoin(aire_suavizado, poblacion_denso, telectricocdmx_denso, vehic; on=:AÑO, makeunique=true)
	#transform!(MAIN, [:CONTAMINANTE, :ALCALDIA] .=> categorical; renamecols=false)
	sort!(MAIN, :AÑO)
end
  ╠═╡ =#

# ╔═╡ c8174734-b36c-4113-90c4-eec8ba5188b3
# CSV.write("procesados.csv", MAIN; bom=true)

# ╔═╡ 57659165-e2d8-47fb-8bb3-3f2b844cb3fa
md"## Población"

# ╔═╡ Cell order:
# ╠═0314c9ba-8e88-11ec-0b13-ff6d7516e4b5
# ╠═84782bc2-92c5-4662-95f9-09c3312b6d4d
# ╟─58cdca1b-0ae8-409c-ade4-272fc741b392
# ╟─b21f5c87-2700-4914-8ab4-8c1190e21148
# ╟─222349ef-d993-4af0-9975-f94acae44b88
# ╠═893be8c3-90e5-4498-bab5-2b500045f989
# ╠═cfc25d9d-35ae-48f8-83e5-bce6c415eea6
# ╠═0ef2d15e-14ae-47cf-b73e-2565b621f670
# ╠═32d95ebc-ae09-4cec-ab64-2e495032a008
# ╟─d930feb1-d7dc-41df-8d53-478e2d1c5a71
# ╠═0ff0079e-24eb-4edd-865c-8b9fd91215ea
# ╠═5110c22b-cebe-460f-9a0a-b65f2643c866
# ╟─bf354eac-7ed5-4713-828e-305de45f502a
# ╠═c28ac129-afd3-4071-b1ac-02d8b1c7565a
# ╠═04a38ace-8c57-401a-8759-43d9d77175d0
# ╠═a27b24b4-e9a5-419e-bcea-25b8cfbffc66
# ╠═d35505d9-2fcd-43a7-9b06-b5d26be1e349
# ╠═491f2be2-799f-42a5-8de7-d7fed574d3ce
# ╠═5540ed04-e09d-4d56-8587-5b2afae5abc2
# ╟─5641e8e1-8b86-4f94-b99c-c6cc63bd18be
# ╠═1bbc9870-431f-4511-9292-62dd30ea1383
# ╠═4a071190-96ba-48c3-a262-5487f4c6155d
# ╠═dc4208eb-e975-4c46-9d06-419e20bd1505
# ╠═af5e5030-fbc4-4e7f-8b90-8fe214892290
# ╠═0cb4cac9-d546-4f1f-ac46-9db51b990baa
# ╟─993cb8bf-aba6-4b3e-888c-d3c639967a52
# ╠═438104ae-88bd-45b9-b87d-9a16938d35c1
# ╠═fa1bbb87-ed47-4ac5-b9b5-7cd511bd8349
# ╠═71209b46-ecf2-4721-9e04-f64e76411b64
# ╠═6ddf350b-b847-4867-afd5-d7bdab5fb62d
# ╠═3f725649-ae2c-41c2-899e-b4bcb9ee273b
# ╠═15060361-9c2f-4557-b45c-1e893381d398
# ╠═c8174734-b36c-4113-90c4-eec8ba5188b3
# ╟─57659165-e2d8-47fb-8bb3-3f2b844cb3fa
