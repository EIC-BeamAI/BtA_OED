using BtaUQ
beam=default_beam_moments()
using BeamTracking
lat=bta_lattice(;apertures=false)
for n in (1,64,1024), simd in (false,true)
 function run_particles()
  b=Bunch(repeat(permutedims(beam.mean),n);species=lat.species_ref,p_over_q_ref=lat.p_over_q_ref)
  BeamTracking.track!(b,lat;use_explicit_SIMD=simd,use_cpu_multithreading=false,rf_on=false)
  b.coords.v
 end
 run_particles()
 t=sort([@elapsed run_particles() for _ in 1:3])[2]
 println("particles=",n," simd=",simd," seconds=",t);flush(stdout)
end
