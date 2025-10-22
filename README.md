# Lunar Gravity Field Calculation Project

![Project Cover](cover.png)

A Fortran project for calculating lunar gravity fields using spherical harmonic models for high-precision gravitational field computation.

## Project Overview

This project provides a complete framework for lunar gravity field calculation, including:

- **Gravity Field Model Definition**: Supports spherical harmonic coefficient loading and gravitational acceleration calculation
- **Modular Design**: Uses object-oriented style Fortran module design
- **High-Precision Calculation**: Supports spherical harmonic expansion up to degree 80 and order 80
- **Test Program**: Includes comprehensive test cases to verify calculation accuracy

## File Structure

```
lunarGravity/
├── module_Gravity.f90      # Core gravity field calculation module
├── test_lunar_gravity.f90  # Test program
├── grail.txt              # Lunar gravity field coefficient data file
├── README.md              # Project documentation (English)
├── README-中文.md         # Project documentation (Chinese)
└── .gitignore            # Git ignore file configuration
```

## Core Features

### Gravity Field Object (GravityField)

The gravity field object contains the following attributes:
- Gravity field name and central body
- GM value (gravitational constant × mass)
- Reference radius
- Maximum degree and order
- Spherical harmonic coefficient matrices C and S
- Coefficient loading status

### Main Methods

1. **Initialization** (`initialize_gravity_field`): Create gravity field object
2. **Read Coefficients** (`read_gravity_coefficients`): Load spherical harmonic coefficients from file
3. **Calculate Acceleration** (`compute_gravity_acceleration`): Calculate gravitational acceleration at given position
4. **Display Information** (`get_gravity_field_info`): Output basic gravity field information

## Compilation and Execution

### Compilation Method

Compile the project using a Fortran compiler:

```bash
# Compile module
gfortran -c module_Gravity.f90 -o module_Gravity.o

# Compile test program
gfortran -c test_lunar_gravity.f90 -o test_lunar_gravity.o

# Link to generate executable
gfortran module_Gravity.o test_lunar_gravity.o -o test_lunar_gravity.exe
```

### Run Test

```bash
./test_lunar_gravity.exe
```

## Data File Format

`grail.txt` file format:
- First line: Reference radius, GM value, reference longitude, reference latitude
- Subsequent lines: Degree l, order m, C coefficient, S coefficient, C uncertainty, S uncertainty

## Technical Features

1. **Spherical Harmonic Calculation**: Uses fully normalized associated Legendre polynomials
2. **Numerical Stability**: Employs recursive algorithms to ensure computational stability
3. **Modular Design**: Easy to extend for supporting gravity fields of other celestial bodies
4. **High Precision**: Supports high-order spherical harmonic expansion, meeting scientific research requirements

## Application Scenarios

- Lunar orbital dynamics analysis
- Lunar probe orbit design
- Lunar gravity field scientific research
- Celestial mechanics teaching demonstrations

## Requirements

- Fortran compiler (gfortran, ifort, etc.)
- Basic linear algebra operation libraries

## License

This project uses an open source license. Please refer to the LICENSE file for specific information.

## Contributing

Welcome to submit Issues and Pull Requests to improve the project.

## Contact

Please contact the project maintainer if you have any questions.
