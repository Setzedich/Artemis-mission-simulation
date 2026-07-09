# Artemis Mission Simulation

An interactive space mission simulator built with **Python**, **Streamlit**, **NumPy**, **Pandas**, and **Plotly**.

This project models a simplified **Artemis II free-return lunar trajectory** using **parametric cubic splines** and numerical methods. It provides an interactive environment to visualize a spacecraft trajectory between Earth and the Moon, analyze real-time telemetry, and export mission data.

---

## Demo

<p align="center">
  <img src="videoartemis.gif" alt="Artemis Mission Simulation Demo" width="850"/>
</p>

---

## Overview

**Artemis Mission Simulation** is an educational and scientific visualization project focused on trajectory modeling, numerical methods, and interactive data visualization.

The simulator allows users to modify orbital parameters, observe the spacecraft path in 3D and 2D, track the current mission status, and review telemetry values such as velocity, acceleration, apparent G-force, distance traveled, and current position.

The trajectory is generated using cubic spline interpolation across the X, Y, and Z axes. The spline system is solved using Gaussian elimination with partial pivoting to create a smooth and continuous spacecraft path.

---

## Features

- Interactive 3D visualization of Earth, Moon, and spacecraft trajectory
- 2D trajectory analysis in the X-Y and X-Z planes
- Real-time mission telemetry dashboard
- Adjustable orbital inclination
- Adjustable trajectory crossing point
- Mission time slider to track spacecraft position
- Current vehicle coordinates in kilometers
- Distance-to-Moon calculation
- Navigation status alerts
- CSV export of telemetry data
- Internal mathematical engine for cubic spline calculations
- Gaussian elimination with partial pivoting

---

## Technologies Used

| Technology | Purpose |
|---|---|
| Python | Core programming language |
| Streamlit | Interactive web application |
| NumPy | Numerical calculations |
| Pandas | Telemetry data handling |
| Plotly | 2D and 3D visualizations |

---

## Mathematical Model

The spacecraft trajectory is modeled using **parametric cubic splines**.

Each coordinate axis is calculated independently:

```text
x(t), y(t), z(t)
```

The simulator builds a tridiagonal system of equations to ensure smoothness and continuity between trajectory segments.

The numerical process includes:

1. Construction of the spline coefficient matrix
2. Construction of the result vectors for each axis
3. Gaussian elimination
4. Partial pivoting
5. Back substitution
6. Cubic spline coefficient calculation
7. Spline evaluation over the mission timeline

This method generates a continuous and smooth approximation of the spacecraft trajectory.

---

## Telemetry

The application calculates and displays mission data such as:

- Relative velocity
- Acceleration
- Apparent G-force
- Total distance traveled
- Current X coordinate
- Current Y coordinate
- Current Z coordinate
- Distance to the Moon
- Navigation status

Telemetry data can also be exported as a `.csv` file for further analysis.

---

## Project Structure

```text
Artemis-Mission-Simulation/
│
├── app_artemis.py        # Main Streamlit application
├── Splines.ipynb         # Notebook related to spline calculations
├── README.md             # Project documentation
└── requirements.txt      # Project dependencies
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/your-username/Artemis-Mission-Simulation.git
```

Enter the project folder:

```bash
cd Artemis-Mission-Simulation
```

Install the required dependencies:

```bash
pip install -r requirements.txt
```

If you do not have a `requirements.txt` file, install the dependencies manually:

```bash
pip install numpy pandas plotly streamlit
```

---

## Usage

Run the Streamlit application:

```bash
streamlit run app_artemis.py
```

You can also run it using:

```bash
python -m streamlit run app_artemis.py
```

After running the command, the application will open in your browser.

---

## Requirements

Create a `requirements.txt` file with the following dependencies:

```txt
numpy
pandas
plotly
streamlit
```

---

## Application Preview

The simulator includes:

- A main 3D orbital viewer
- Earth and Moon rendered as 3D objects
- Artemis II trajectory path
- Orion spacecraft position marker
- 2D orbital analysis graphs
- Real-time telemetry metrics
- Mission status alerts
- Downloadable telemetry table
- Internal mathematical model viewer

---

## Educational Purpose

This project was developed to apply and demonstrate concepts from:

- Numerical methods
- Cubic spline interpolation
- Linear systems
- Scientific computing
- Space trajectory visualization
- Interactive data applications

It combines mathematical modeling with modern visualization tools to create an intuitive simulation environment.

---

## Author

**Marlon Molina**  
Cybernetics Engineering in Computer Systems Student

---

## License

This project is licensed under the MIT License.
