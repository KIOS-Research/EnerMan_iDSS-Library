%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 1: A library for hydraulic calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 1.1: constants
%%%%%%%%%%%%%%%%%%%%%%%%%%

% g (m/s^2)
constant(g, 9.81).
constant(pi, 3.1416).

% rho (Kg/m^3)
constant(rho, 1000).
% mu Kg/(m*s)
constant(mu, 0.000547).

fluid(fluid(Density, Viscocity)) :-
    constant(rho, Density),
    constant(mu, Viscocity).

%%%%%%%%%%%%%%%%%%%%%%%%%%
% section 1.2: poisson factors
%%%%%%%%%%%%%%%%%%%%%%%%%%

%poisson(rubber, 0.48 - ~0.5)

poisson(lead, 0.431).

poisson(clay, 0.41).
poisson(sandy_clay, 0.37).
poisson(z_nickel, 0.36).
poisson(phosphor_bronze, 0.359).
poisson(brass_cast, 0.357).
poisson(copper, 0.355).
poisson(aluminum_6061_T6, 0.35).
poisson(magnesium, 0.35).
poisson(bronze, 0.34).
poisson(polystyrene, 0.34).
poisson(aluminum, 0.334).
poisson(brass_70_30, 0.331).
poisson(zinc, 0.331).
poisson(ice, 0.33).
poisson(nickel_silver, 0.322).
poisson(aluminum_2024_T4, 0.32).
poisson(titanium_99_0_Ti, 0.32).
poisson(monel_metal, 0.315).
poisson(sandy_loam, 0.31).
poisson(molybdenum, 0.307).
poisson(stainless_steel_18_8, 0.305).
poisson(steel_mild, 0.303).
poisson(steel_high_carbon, 0.295).
poisson(nickel_steel, 0.291).
poisson(sand, 0.29).
poisson(steel_cold_rolled, 0.287).
poisson(beryllium_copper, 0.285).
poisson(magnesium_alloy, 0.281).
poisson(wrought_iron, 0.278).
poisson(iron_malleable, 0.271).
%poisson(inconel, 0.27 - 0.38).
poisson(steel_cast, 0.265).
poisson(iron_ductile, 0.28).
poisson(glass_soda, 0.22).
%poisson(iron_cast, 0.22 - 0.30).
poisson(iron_cast_gray, 0.211).
%poisson(glass_float, 0.2 - 0.27).
%poisson(granite, 0.2 - 0.3).
%poisson(limestone, 0.2 - 0.3).
%poisson(marble, 0.2 - 0.3).
%poisson(concrete, 0.1 - 0.2).
poisson(cork, 0).
poisson(upvc, 0.34).

% https://www.neonalloys.com/304-stainless-steel
poisson(stainless_steel_304SS, 0.270).
% https://super-metals.com/wp-content/uploads/2015/03/SS-316.pdf
poisson(stainless_steel_316SS, 0.270).
% https://www.azom.com/properties.aspx?ArticleID=712
poisson(titanium, 0.36).

poisson(Specs, Poisson_Factor) :-
    material(Specs, Material),
    poisson(Material, Poisson_Factor).

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 1.3: component specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%

%    class model: specs -> pipe_specs(Diameter (m), Length (m), Height1 (m), Height2 (m), Material)
%    class model: specs -> heat_exchanger_specs(Length (m), Width (m), Distance (m), Material)

component_length(pipe_specs(_, Length, _, _, _), Length).
component_length(heat_exchanger_specs(Length, _, _, _), Length).

plate_width(heat_exchanger_specs(_, Width, _, _), Width).
plate_distance(heat_exchanger_specs(_, _, Distance, _), Distance).

material(pipe_specs(_, _, _, _, Material), Material).
material(heat_exchanger_specs(_, _, _, Material), Material).

height_difference(pipe_specs(_, _, Height1, Height2, _), DH) :-
    DH is Height2 - Height1.

height_difference(heat_exchanger_specs(_, _, _, _), 0).

% hydraulic_diameter: meter
hydraulic_diameter(pipe_specs(Diameter, _, _, _, _), Diameter).

% hydraulic_diameter: meter
% Dh = hydraulic diameter = [4*l*dplate] / [2*(l+dplate)] (m)
% l = width of the plates (m)
% dplate = gap in between 2 plates (m)
hydraulic_diameter(HeatExchangerSpecs, Diameter) :-
    plate_width(HeatExchangerSpecs, Width),
    plate_distance(HeatExchangerSpecs, Distance),
    Diameter is (2 * Width * Distance) / (Width + Distance).

% area: m^2
area(pipe_specs(Diameter, _, _, _, _), Area) :-
    constant(pi, Pi),
    Area is (Pi * Diameter * Diameter) / 4.

% area: m^2
area(heat_exchanger_specs(L, Width, D, _), Area) :-
    hydraulic_diameter(heat_exchanger_specs(L, Width, D, _), Diameter),
    Area is Diameter * Width.

% this is for computing the friction_head_loss
% (there is a different coefficient between pipes and heat echangers)
component_coef(pipe_specs(_, _, _, _, _), 1).
component_coef(heat_exchanger_specs(_, _, _, _), 4).

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 1.4: friction factor
%%%%%%%%%%%%%%%%%%%%%%%%%%

reynolds_number(Flow, fluid(Density, Viscosity), Specs, Reynolds) :-
    hydraulic_diameter(Specs, Diameter),
    Reynolds is (Density * Flow * Diameter) / Viscosity.

component_roughness(Specs, Roughness) :-
    hydraulic_diameter(Specs, Diameter),
    poisson(Specs, Poisson_Factor),
    Roughness is Poisson_Factor / Diameter.

friction_factor(Flow, fluid(Density, Viscosity), Specs, Friction_Factor) :-
    reynolds_number(Flow, fluid(Density, Viscosity), Specs, Reynolds),
    (
            Reynolds =< 2300
        ->
            Friction_Factor is 64 / Reynolds
        ;
            (
                component_roughness(Specs, Roughness),
                Friction_Factor is 0.25 / (log10((Roughness / 3.7) + (5.74 / (Reynolds^0.9)))^2)
            )
    ).

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 1.5: friction loss and head loss
%%%%%%%%%%%%%%%%%%%%%%%%%%

% velocity: m/s
velocity(Flow, Specs, Velocity) :-
    area(Specs, Area),
    Velocity is Flow / Area.

% flow: m * m * m /s
flow(Velocity, Specs, Flow) :-
    area(Specs, Area),
    Flow is Velocity * Area.

% friction_head_loss: m
friction_head_loss(Flow, fluid(Density, Viscosity), Specs, Friction_Head_Loss) :-
    friction_factor(Flow, fluid(Density, Viscosity), Specs, Friction_Factor),
    component_length(Specs, Length),
    hydraulic_diameter(Specs, Diameter),
    velocity(Flow, Specs, Velocity),
    constant(g, G),
    component_coef(Specs, Component_Coef),
    Friction_Head_Loss is Component_Coef * (Friction_Factor * Velocity * Velocity * Length) / (2 * Diameter * G).

%%%%%%%%%%%%%%%%%%%%%%%%%%
% pressure drop in a heat HeatExchanger
% ΔP = 4 f ρ u u L / (2 Dh)
%%%%%%%%%%%%%%%%%%%%%%%%%%

%friction_head_loss(Flow, fluid(Density, Viscosity), HeatExchangerSpecs, Friction_Head_Loss) :-
%    friction_factor(Flow, fluid(Density, Viscosity), HeatExchangerSpecs, Friction_Factor),
%    component_length(HeatExchangerSpecs, Length),
%    hydraulic_diameter(HeatExchangerSpecs, Diameter),
%    constant(g, G),
%    Friction_Head_Loss is 4 * (Friction_Factor * Density * Velocity * Velocity * Length) / (2 * Diameter * G).

%%%%%%%%%%%%%%%%%%%%%%%%%%
% pressure drop in a pipe
%
%    P1 + (1/2 * ρ * v1 * v1) + (ρ * g * h1) = P2 + (1/2 * ρ * v2 * v2) + (ρ * g * h2) + friction_head_loss
% => P1 - P2 = (1/2 * ρ * v2 * v2) + (ρ * g * h2) + friction_head_loss - (1/2 * ρ * v1 * v1) - (ρ * g * h1)
% => ΔP = 1/2 * ρ * (v2 * v2 - v1 * v1) + ρ * g * (h2 - h1) + friction_head_loss
%%%%%%%%%%%%%%%%%%%%%%%%%%
% head_loss (m): Flow (m^3/s), Density (Kg/m^3), Viscocity ((Kg/m*s))
head_loss(Flow, fluid(Density, Viscosity), Specs, Head_Loss) :-
    friction_head_loss(Flow, fluid(Density, Viscosity), Specs, Friction_Head_Loss),
    height_difference(Specs, DiffH),
    Head_Loss is DiffH + Friction_Head_Loss.
%    Pressure_Drop is Rho * G * (Height_2 - Height_1) + (Rho * G * (V2^2 - V1^2) / 2) + Friction_Head_Loss.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 1.6: pump operation
%%%%%%%%%%%%%%%%%%%%%%%%%%

% pump_efficiency (%): Flow (m^3/h), Density (Kg/m^3), Head (m), Consuption (KW)
% pump_efficiency(Flow, Density, Head, Consumption, Efficiency) :-
%    constant(g, G),
%    Efficiency is ((Density * G * (Flow / 3600) * Head) / (Consumption * 1000)).

% pump_efficiency (%): Flow (m^3/s), Density (Kg/m^3), Head (m), Consuption (W)
pump_efficiency(Flow, Density, Head, Consumption, Efficiency) :-
    constant(g, G),
    Efficiency is ((Density * G * Flow * Head) / Consumption).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 2: A hydraulic ontology for the iDSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 2.1; class models
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%
% Section 2.1.1: pipe class models
%%%%%

class_model(pipe, Pipe, head_loss, [Flow, fluid(Density, Viscosity)], Head_Loss) :-
    compute_model(Pipe, pipe_specs, [], PipeSpecs),
    head_loss(Flow, fluid(Density, Viscosity), PipeSpecs, Head_Loss).

%%%%%
% Section 2.1.2: tank class models
%%%%%

class_model(tank, Tank, pressure, [], Pressure) :-
    compute_model(Tank, tank_specs, [], tank_specs(_, Pressure)).

%%%%%
% Section 2.1.3: heat_exchanger class models
%%%%%

class_model(heat_exchanger, HeatExchanger, head_loss, [Flow, fluid(Density, Viscosity)], Head_Loss) :-
    compute_model(HeatExchanger, heat_exchanger_specs, [], Specs),
    head_loss(Flow, fluid(Density, Viscosity), Specs, Head_Loss).
%%%%%
% Section 2.1.4: pump class models
%%%%%

class_model(pump, Pump, energy_efficiency, [Flow, Density], Efficiency) :-
    compute_model(Pump, performance, [Flow], Head),
    compute_model(Pump, consumption, [Flow], Consumption),
    pump_efficiency(Flow, Density, Head, Consumption, Efficiency).


%%%%%
% Section 2.1.4: pump class models
%%%%%

class_model(pump_system, PumpSystem, energy_efficiency, [Flow, Density], Efficiency) :-
    compute_model(PumpSystem, performance, [Flow], Head),
    compute_model(PumpSystem, consumption, [Flow], Consumption),
    pump_efficiency(Flow, Density, Head, Consumption, Efficiency).

%%%%%
% Section 2.1.5: degrease_tank class models
%%%%%

class_model(degrease_tank, Degrease_Tank, recirculation_system_head_loss, [Flow, fluid(Density, Viscosity)], Recirc_head_loss) :-
    compute_model(Degrease_Tank, p1, head_loss, [Flow, fluid(Density, Viscosity)], PD1),
    compute_model(Degrease_Tank, p2, head_loss, [Flow, fluid(Density, Viscosity)], PD2),
    compute_model(Degrease_Tank, p3, head_loss, [Flow, fluid(Density, Viscosity)], PD3),
    compute_model(Degrease_Tank, heat_exchanger, head_loss, [Flow, fluid(Density, Viscosity)], PDHE),
    Recirc_head_loss is (PD1 + PD2 + PD3 + PDHE).


class_model(degrease_tank, Degrease_Tank, feasible_recirculation_pump_system, [Flow, fluid(Density, Viscosity)], Head) :-
    compute_model(Degrease_Tank, recirculation_system_pressure_loss, [Flow, fluid(Density, Viscosity)], Pressure_Loss),
    compute_model(Degrease_Tank, recirculation_pump_system, performance, [Flow], PumpHead),
    (
            PumpHead >= Pressure_Loss
        ->  Head is PumpHead
        ;   undef_value(Head) % is undef %undef(Head)
    ).

class_model(degrease_tank, Degrease_Tank, feasible_recirculation_pump_system_distribution, [Distribution, fluid(Density, Viscosity)], Head) :-
    distribution(Degrease_Tank, feasible_recirculation_pump_system, [fluid(Density, Viscosity)], Distribution, Heads),
    avg(Heads, Head).

class_model(degrease_tank, Degrease_Tank, spray_system_head_loss, [Flow, fluid(Density, Viscosity)], Spray_head_loss) :-
    compute_model(Degrease_Tank, p4, head_loss, [Flow, fluid(Density, Viscosity)], PD4),
    compute_model(Degrease_Tank, p5, head_loss, [Flow, fluid(Density, Viscosity)], PD5),
    compute_model(Degrease_Tank, p7, head_loss, [Flow, fluid(Density, Viscosity)], PD7),
    Spray_head_loss is (PD4 + PD5 + PD7).

class_model(degrease_tank, Degrease_Tank, feasible_spray_pump_system, [Flow, fluid(Density, Viscosity)], Head) :-
    compute_model(Degrease_Tank, spray_system_head_loss, [Flow, fluid(Density, Viscosity)], Pressure_Loss),
    compute_model(Degrease_Tank, spray_system, performance, [Flow], PumpHead),
    (
            PumpHead >= Pressure_Loss
        ->  Head is PumpHead
        ;   undef_value(Head) % is undef %undef(Head)
    ).

class_model(degrease_tank, Degrease_Tank, efficiency_spray_pump_system, [Flow, fluid(Density, Viscosity)], Efficiency) :-
    compute_model(Degrease_Tank, spray_system_head_loss, [Flow, fluid(Density, Viscosity)], Pressure_Loss),
    compute_model(Degrease_Tank, spray_system, performance, [Flow], PumpHead),
    (
            PumpHead >= Pressure_Loss
        ->  compute_model(Degrease_Tank, spray_system, energy_efficiency, [Flow, Density], Efficiency)
        ;   undef_value(Efficiency) % is undef %undef(Efficiency)
    ).

class_model(degrease_tank, Degrease_Tank, energy_efficiency, [Flow, fluid(Density, Viscosity)], Efficiency) :-
  class_model(degrease_tank, Degrease_Tank, efficiency_spray_pump_system, [Flow, fluid(Density, Viscosity)], Efficiency).

%%%%%
% Section 2.1.6: testbed class models
%%%%%

class_model(wdn, Testbed, energy_efficiency_duty_point_distribution, Input, Efficiency) :-
    compute_model(Testbed, b, energy_efficiency_duty_point_distribution, Input, Efficiency).


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section 2.2: model implementations
%%%%%%%%%%%%%%%%%%%%%%%%%%

model(Pump, Polynomial, [Flow], Head) :-
    class(Pump, pump),
    polynomial(Polynomial, Flow, Head).

%%%%%
% Section 2.2.1: pipe models
%%%%%

model(Pipe, pipe_specs(Diameter, Length, Height1, Height2, Material), [], pipe_specs(Diameter, Length, Height1, Height2, Material)) :-
    class(Pipe, pipe).

%%%%%
% Section 2.2.2: tank models
%%%%%

model(Tank, tank_specs(Height, Level), [], tank_specs(Height, Level)) :-
    class(Tank, tank).

%%%%%
% Section 2.2.3: heat_exchanger models
%%%%%

model(HeatExchanger, heat_exchanger_specs(Length, Width, Distance, Material), [], heat_exchanger_specs(Length, Width, Distance, Material)) :-
    class(HeatExchanger, heat_exchanger).

%%%%%
% Section 2.2.4: pump system models
%%%%%

model(Pump_system, single_pump_performance, [Flow], Head) :-
    compute_model(Pump_system, pump_1, performance, [Flow], Head).

model(Pump_system, single_pump_consumption, [Flow], Head) :-
    compute_model(Pump_system, pump_1, consumption, [Flow], Head).




model(Pump_system, duty_standby_performance, [Flow], Head) :-
    compute_model(Pump_system, pump_1, performance, [Flow], Head),
    compute_model(Pump_system, pump_2, performance, [Flow], Head).

model(Pump_system, duty_standby_consumption, [Flow], Head) :-
    compute_model(Pump_system, pump_1, consumption, [Flow], Head),
    compute_model(Pump_system, pump_2, consumption, [Flow], Head).




model(Pump_system, duty_assist_standby_performance, [Flow], Head) :-
    %modelImpl(Pump_system, pump_1, performance, FunctionImpl1),
    %modelImpl(Pump_system, pump_2, performance, FunctionImpl2),
    %addOnX(FunctionImpl1, FunctionImpl2, Flow, Head).
    modelImpl(Pump_system, pump_1, performance, coefficients([A, B, C])),
    %modelImpl(Pump_system, pump_2, performance, coefficients([A, B, C])),
    AA is A / 4,
    BB is B / 2,
    polynomial(coefficients([AA, BB, C]), Flow, Head).


%model(Pump_system, duty_assist_standby_performance, [_], undef) :-
%    modelImpl(Pump_system, pump_1, performance, coefficients([A, B, C])),
%    modelImpl(Pump_system, pump_2, performance, coefficients([AA, BB, CC])),
%    \+ (AA is A, BB is B, CC is C).

model(Pump_system, duty_assist_standby_consumption, [Flow], Consumption) :-
    WaterflowAtPump is Flow / 2,
    compute_model(Pump_system, pump_1, consumption, [WaterflowAtPump], ConsumptionAtPump),
    %compute_model(Pump_system, pump_2, consumption, [WaterflowAtPump], ConsumptionAtPump),
    Consumption is 2 * ConsumptionAtPump. %1 + Consumption2.

model(Pump_system, duty_assist_standby_energy_efficiency, [Flow, Density], Efficiency) :-
    compute_class_model(Pump_system, energy_efficiency, [Flow, Density], Eff1),
    compute_model(Pump_system, pump_1, energy_efficiency, [Flow, Density], Eff2),
    foldr(max_elements, undef, [Eff1, Eff2], Efficiency).
