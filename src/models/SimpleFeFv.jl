struct SimpleFeFv{NP, NS, E} <: AbstractMechanicalModel{NP, NS}
  equilibrium::E
end

function SimpleFeFv(equilibrium)
  NP = num_properties(equilibrium)
  NS = 0
  return SimpleFeFv{NP, NS, typeof(equilibrium)}(equilibrium)
end

