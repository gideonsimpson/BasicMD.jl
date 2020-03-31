abstract type AbstractSampler end
abstract type ZerothOrderSampler <: AbstractSampler
abstract type FirstOrderSampler <: AbstractSampler
abstract type SecondOrderSampler <: AbstractSampler

abstract type AbstractSamplerState end
abstract type ZerothOrderSamplerState <: AbstractSamplerState
abstract type FirstOrderSamplerState <: AbstractSamplerState
abstract type SecondOrderSamplerState <: AbstractSamplerState
