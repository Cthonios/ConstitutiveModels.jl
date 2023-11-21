function strip_units(qoi)
  return ustrip(qoi), unit(qoi)
end

function strip_units(qoi::A) where {A <: AbstractArray}
  return ustrip.(qoi), unit.(qoi)
end

function add_units(qoi, unit_to_add)
  qoi * unit_to_add
end

function add_units(qoi::A1, unit_to_add::A2) where {A1 <: AbstractArray, A2 <: AbstractArray}
  qoi .* unit_to_add
end