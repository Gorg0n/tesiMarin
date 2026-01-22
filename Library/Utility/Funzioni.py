import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
from Library.AuxiliaryComponents.Constants import *


def compute_h2_demand(yearly_h2_demand, year_length, time):
    '''
    :yearly_h2_demand:float -> yearly hydrogen demand [ton/year]
    :year_length:float -> length of the year [h/year]
    :time:string -> time of discharge (hourly, half-day, daily, weekly)
    OUTPUT:
        h2_demand:list -> list of hydrogen demand for each time step (kg/h, kg/day, kg/week)
    '''
    hourly_h2_demand = yearly_h2_demand*1e3/year_length # [kg/h]

    h2_demand = []
    for i in range(year_length): 
        if time == 'hourly':
            h2_demand.append(hourly_h2_demand) 
        elif time == 'half-day':
            if i % 12 == 0:
                h2_demand.append(hourly_h2_demand*12)
            else:
                h2_demand.append(0)
        elif time == 'daily':
            if i % 24 == 0:
                h2_demand.append(hourly_h2_demand*24)
            else:
                h2_demand.append(0)
        elif time == 'weekly':
            if i % 168 == 0:
                h2_demand.append(hourly_h2_demand*168)
            else:
                h2_demand.append(0)

    return h2_demand


def compute_hydrogen_density(P, T): 
    '''
    :input: P:float -> pressione [bar]
    :input: T:float -> temperatura [°C]
    :input: z:float -> fattore di compressibilità (default 1.0 → gas ideale)
    reference: Standardized Equation for Hydrogen Gas Densities for Fuel Consumption Applications1 - Lemmon
    '''
    P = P * conv_bar_pascal  # conversione bar → Pascal
    T = T + conv_celsius_kelvin  # conversione °C → K
    coeffs = [
    {"i": 2, "j": 1, "v_ij": 3.6719e-2, "n_ij": -1.23},
    {"i": 2, "j": 2, "v_ij": -3.9839e-2, "n_ij": -2.22},
    {"i": 3, "j": 1, "v_ij": -1.4722e-3, "n_ij": -2.68},
    {"i": 3, "j": 2, "v_ij": 2.4083e-3, "n_ij": -3.1},
    {"i": 4, "j": 1, "v_ij": 0.65994e-5, "n_ij": -2.7},
    {"i": 4, "j": 2, "v_ij": -0.15469e-4, "n_ij": -4.3},
    {"i": 5, "j": 2, "v_ij": -0.13383e-6, "n_ij": -3.3},
    {"i": 6, "j": 1, "v_ij": 0.15608e-8, "n_ij": -4.1},
    ]
    z = 1 # fattore di compressibilità (default 1.0 → gas ideale)
    for term in coeffs:
        i = term["i"]
        j = term["j"]
        v_ij = term["v_ij"]
        n_ij = term["n_ij"]
        z += v_ij * (P*1e-6) ** (i - 1) * (T / 100.0) ** n_ij
    rho = (P * massa_molare_h2) / (z * r_gas * T)  # densità [kg/m³]

    return rho 

def plot_multiple_curves(y_set_values, x_values, labels, xlabel, ylabel, linestyle_set=None): 
    '''
    y_set_values:list of lists -> ogni lista contiene i valori della specifica curva
    x_values:list -> valori dell'asse x 
    labels:list of str -> etichette per ogni curva
    '''
    plt.figure(figsize=(10, 6))
    if linestyle_set is None:
        for y, label in zip(y_set_values, labels):
            plt.plot(x_values, y, label=label)
    else: 
        for y, label, line_style in zip(y_set_values, labels, linestyle_set):
            plt.plot(x_values, y, label=label, linestyle=line_style)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=16)

def scatter_plot_multiple_curves(y_set_values, x_values, labels, xlabel, ylabel): 
    '''
    y_set_values:list of lists -> ogni lista contiene i valori della specifica curva
    x_values:list -> valori dell'asse x 
    labels:list of str -> etichette per ogni curva
    '''
    plt.figure(figsize=(10, 6))
    for y, label in zip(y_set_values, labels):
        plt.scatter(x_values, y, label=label)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=16)

def plot_single_curve(y_values, x_values, xlabel, ylabel, linestyle=None): 
    '''
    y_values:list -> valori dell'asse y
    x_values:list -> valori dell'asse x 
    '''
    plt.figure(figsize=(10, 6))
    if linestyle is None:
        plt.plot(x_values, y_values)
    else: 
        plt.plot(x_values, y_values, linestyle=linestyle)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.grid(True)
    plt.tick_params(axis='both', which='major', labelsize=16)

