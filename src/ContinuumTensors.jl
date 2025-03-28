abstract type AbstractContinuumTensor{T} end

function Base.show(io::IO, A::AbstractContinuumTensor)
  println(io, "$(typeof(A)):")
  display(A.val)
end

@inline Base.getindex(A::AbstractContinuumTensor, i::Int, j::Int) = A.val[i, j]

struct TwoPointTensor{T} <: AbstractContinuumTensor{T}
  val::Tensor{2, 3, T, 9}
end

@inline function TwoPointTensor(vals...)
  return TwoPointTensor(Tensor{2, 3, eltype(vals), 9}((vals)))
end

@inline function TwoPointTensor(A::TwoPointTensor)
  return TwoPointTensor(A.val)
end

@inline function TwoPointTensor(A::SMatrix{3, 3, T, 9}) where T
  return TwoPointTensor(A...)
end

@inline Base.:+(A::TwoPointTensor, B::TwoPointTensor) = TwoPointTensor(A.val + B.val)
@inline Base.:-(A::TwoPointTensor) = TwoPointTensor(-A.val)
@inline Base.:-(A::TwoPointTensor, B::TwoPointTensor) = TwoPointTensor(A.val - B.val)
@inline Base.:*(c::T, A::TwoPointTensor{T}) where T <: Number = TwoPointTensor{T}(c * A.val)
@inline Tensors.det(A::TwoPointTensor) = det(A.val)
@inline Tensors.one(A::TwoPointTensor) = TwoPointTensor(one(A.val))
@inline Tensors.tovoigt(type, A::TwoPointTensor) = tovoigt(type, A.val)
@inline Tensors.tr(A::TwoPointTensor) = tr(A.val)

struct TwoPointTensorTranspose{T} <: AbstractContinuumTensor{T}
  val::Tensor{2, 3, T, 9}
end

@inline Base.:+(A::TwoPointTensorTranspose, B::TwoPointTensorTranspose) = TwoPointTensorTranspose(A.val + B.val)
@inline Base.:-(A::TwoPointTensorTranspose) = TwoPointTensorTranspose(-A.val)
@inline Base.:-(A::TwoPointTensorTranspose, B::TwoPointTensorTranspose) = TwoPointTensorTranspose(A.val - B.val)
@inline Base.:*(c::T, A::TwoPointTensorTranspose{U}) where {T <: Number, U} = TwoPointTensorTranspose{U}(c * A.val)
@inline Tensors.det(A::TwoPointTensorTranspose) = det(A.val)
@inline Tensors.one(A::TwoPointTensorTranspose) = TwoPointTensorTranspose(one(A.val))
@inline Tensors.tr(A::TwoPointTensorTranspose) = tr(A.val)

# TODO check on these, especially the inverse
@inline Tensors.adjoint(A::TwoPointTensor) = TwoPointTensorTranspose(A.val')
@inline Tensors.adjoint(A::TwoPointTensorTranspose) = TwoPointTensor(A.val')
@inline Tensors.inv(A::TwoPointTensor) = TwoPointTensorTranspose(inv(A.val))
@inline Tensors.inv(A::TwoPointTensorTranspose) = TwoPointTensorTranspose(inv(A.val))
@inline Tensors.transpose(A::TwoPointTensor) = TwoPointTensorTranspose(A.val')
@inline Tensors.transpose(A::TwoPointTensorTranspose) = TwoPointTensor(A.val')

struct MaterialTensor{T} <: AbstractContinuumTensor{T}
  val::SymmetricTensor{2, 3, T, 6}
end

@inline Base.:+(A::MaterialTensor, B::MaterialTensor) = MaterialTensor(A.val + B.val)
@inline Base.:-(A::MaterialTensor) = MaterialTensor(-A.val)
@inline Base.:-(A::MaterialTensor, B::MaterialTensor) = MaterialTensor(A.val - B.val)
@inline Base.:*(c::T, A::MaterialTensor{U}) where {T <: Number, U} = MaterialTensor{U}(c * A.val)
@inline Tensors.adjoint(A::MaterialTensor) = MaterialTensor(A.val)
@inline Tensors.dev(A::MaterialTensor) = MaterialTensor(dev(A.val))
@inline Tensors.dcontract(A::MaterialTensor, B::MaterialTensor) = dcontract(A.val, B.val)
@inline Tensors.norm(A::MaterialTensor) = norm(A.val)
@inline Tensors.one(A::MaterialTensor) = MaterialTensor(one(A.val))
@inline Tensors.tovoigt(type, A::MaterialTensor) = tovoigt(type, A.val)
@inline Tensors.tr(A::MaterialTensor) = tr(A.val)
@inline Tensors.transpose(A::MaterialTensor) = MaterialTensor(A.val)
@inline Tensors.tdot(A::TwoPointTensor) = MaterialTensor(Tensors.tdot(A.val))

struct SpatialTensor{T} <: AbstractContinuumTensor{T}
  val::SymmetricTensor{2, 3, T, 6}
end

function SpatialTensor(vals...)
  return SpatialTensor(SymmetricTensor{2, 3, eltype(vals), 6}((vals)))
end

@inline Base.:+(A::SpatialTensor, B::SpatialTensor) = SpatialTensor(A.val + B.val)
@inline Base.:-(A::SpatialTensor) = SpatialTensor(-A.val)
@inline Base.:-(A::SpatialTensor, B::SpatialTensor) = SpatialTensor(A.val - B.val)
@inline Base.:*(c::T, A::SpatialTensor{U}) where {T <: Number, U} = SpatialTensor{U}(c * A.val)
@inline Tensors.adjoint(A::SpatialTensor) = SpatialTensor(A.val)
@inline Tensors.dcontract(A::SpatialTensor, B::SpatialTensor) = dcontract(A.val, B.val)
@inline Tensors.dev(A::SpatialTensor) = SpatialTensor(dev(A.val))
@inline Tensors.dott(A::TwoPointTensor) = SpatialTensor(Tensors.dott(A.val))
@inline Tensors.norm(A::TwoPointTensor) = norm(A.val)
@inline Tensors.one(A::SpatialTensor) = SpatialTensor(one(A.val))
@inline Tensors.tovoigt(type, A::SpatialTensor) = tovoigt(type, A.val)
@inline Tensors.tr(A::SpatialTensor) = tr(A.val)
@inline Tensors.transpose(A::SpatialTensor) = SpatialTensor(A.val)

# TODO check this one
@inline Base.:*(A::SpatialTensor, B::TwoPointTensor) = TwoPointTensor(dot(A.val, B.val))

@inline Base.:*(A::TwoPointTensor, B::TwoPointTensorTranspose) = SpatialTensor(symmetric(dot(A.val, B.val)))
