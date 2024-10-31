import scipy as sp
import numpy as np
import pylab as plt
import random
from matplotlib.animation import FuncAnimation
from scipy.optimize import rosen_hess, rosen_der
from scipy.optimize import rosen
import matplotlib.cm as cm
#
# Define the 2-D Rosenbrock function
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def brute_force(step):
    x = np.linspace(-2, 3, step)
    X, Y = np.meshgrid(x, x)
    ros = rosen([X, Y])
    ## Brute force
    print(ros.shape)
    tmp = ros[0,0]
    min = ros[0,0]
    min_point = [0, 0]
    a_i = 0
    a_j = 0
    for i in range(ros.shape[0]):
        for j in range(ros.shape[1]):
            tmp = ros[i, j]
            if tmp < min:
                min = tmp
                a_i = i
                a_j = j
                min_point = [x[i], x[j]]

    real_min = ros[a_i, a_j]

    ax = plt.subplot(111, projection='3d')
    ax.plot_surface(X, Y, ros)
    plt.show()
    return min_point, real_minf


def gradient_descent_rosen(initial_point=(-1, 1), step_size=1e-3, tol=1e-10, max_iter=10000):
    x, y = initial_point  # Start at the given initial point
    points = [(x, y)]  # To store the path of points

    for i in range(max_iter):
        grad = rosen_der([x, y])  # Gradient at current point
        x, y = x - step_size * grad[0], y - step_size * grad[1]
        points.append((x, y))

        # Check if we're close to the minimum
        if np.linalg.norm(grad) < tol:
            break

    return points, (x, y), rosen([x, y])

def newtons_method_rosen(initial_point=(-1, 1), tol=1e-10, max_iter=10000):
    x, y = initial_point  # Start at the given initial point
    points = [(x, y)]  # To store the path of points

    for i in range(max_iter):
        grad = rosen_der([x, y])
        hess = rosen_hess([x, y])

        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            print("Hessian is singular, can't proceed with Newton's Method")
            break

        x, y = x - step[0], y - step[1]
        points.append((x, y))

        # Check if close enough to minimum
        if np.linalg.norm(grad) < tol:
            break

    return points, (x, y), rosen([x, y])

def plot_rosenbrock_paths(all_points, labels):
    # Generate grid for surface plot
    x_vals = np.linspace(-2, 3, 1000)
    X, Y = np.meshgrid(x_vals, x_vals)
    Z = rosen([X, Y])

    plt.figure(figsize=(10, 6))
    plt.contourf(X, Y, Z, levels=5)
    plt.colorbar()
    plt.contour(X, Y, Z, levels=[1, 10, 100, 500, 1000], colors="white", linestyles="--")

    # Plot each path with a unique color
    for points, label in zip(all_points, labels):
        points = np.array(points)
        plt.plot(points[:, 0], points[:, 1], '-', label=label)

    plt.plot(1, 1, 'bo', label="Global Minimum")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

# Define different starting points
initial_points = [(-1, 0), (0, -0.5), (1, 2)]

# Run gradient descent from different starting points
gd_paths = []
for i, initial in enumerate(initial_points):
    points, min_point, min_value = gradient_descent_rosen(initial)
    gd_paths.append(points)
    print(f"Gradient Descent - Start {initial}: Min Point: {min_point}, Min Value: {min_value}")

# Run Newton's Method from different starting points
nm_paths = []
for i, initial in enumerate(initial_points):
    points, min_point, min_value = newtons_method_rosen(initial)
    nm_paths.append(points)
    print(f"Newton's Method - Start {initial}: Min Point: {min_point}, Min Value: {min_value}")

# Plot all paths for gradient descent and Newton's method
# plot_rosenbrock_paths(gd_paths, labels=[f"Gradient Descent Start {pt}" for pt in initial_points])
# plot_rosenbrock_paths(nm_paths, labels=[f"Newton's Method Start {pt}" for pt in initial_points])


def comp_area(points):
    if points.shape[0] != 3:
        raise ValueError("The simplex must have exactly 3 points to form a triangle.")

    # Calculate the distances (side lengths) between the points
    a = np.linalg.norm(points[0] - points[1])
    b = np.linalg.norm(points[1] - points[2])
    c = np.linalg.norm(points[2] - points[0])

    # Calculate the semi-perimeter
    s = (a + b + c) / 2

    # Apply Heron's formula for the area of a triangle
    area = np.sqrt(s * (s - a) * (s - b) * (s - c))

    print(f"Calculated area: {area}")
    return area
def reflect(points, p_vals, ros):
    ## Find min
    p_min = np.min(p_vals)
    min_loc = np.where(p_vals == p_min)[0][0]
    p1 = points[min_loc]
    reflection = [0] * len(points)
    diff_list = []
    for i in range(len(points)):
        ## Find the difference between points the points
        p_i = points[i]
        ## p1 to p_i
        if np.all(p1 == p_i):
            continue
        else:
            # calculate the reflection point
            reflection[i] = p_i ## these dont change
            diff = -1 * (p1 - p_i)
            diff_list.append(diff)
    p_reflect = points[1] + diff_list[0]
    reflection[0] = p_reflect
    r_vals = []
    for p in reflection:
        r_vals.append(ros[p[0], p[1]])
    reflection = np.array(reflection)
    r_vals = np.array(r_vals)
    return reflection, r_vals

def grow(points, p_vals, grow_units, ros):
    ## I want to grow the entire triangle to see if there exist another min just outside that triangle
    ## Find centorid
    centroid = np.sum(points, axis = 0)//len(points) ## This will be pointing to the lattice.
    print(centroid)
    grown_points = []
    for p in points:
        # Calculate the direction vector from the centroid to the point
        direction = p - centroid
        if np.linalg.norm(direction) != 0:
            direction = direction // np.linalg.norm(direction) ## This normalizes the direction vector.
            ## Such that it grows outward of the center point
        grown_point = (p + grow_units * direction).astype(int)
        grown_points.append(grown_point)

    g_vals = []
    for p in grown_points:
        g_vals.append(ros[p[0], p[1]])

    g_vals = np.array(g_vals)
    grown_points = np.array(grown_points)

    return grown_points, g_vals

def shrink(points, p_vals, shrink_units, ros):
    ## I want to shrink the entire triangle to see if there exist another min just outside that triangle
    ## Find centorid
    centroid = np.sum(points, axis = 0)//len(points) ## This will be pointing to the lattice.
    shrink_points = []
    for p in points:
        # Calculate the direction vector from the centroid to the point
        direction = centroid - p
        if np.linalg.norm(direction) != 0:
            direction = direction // np.linalg.norm(direction) ## This normalizes the direction vector.
            ## Such that it grows outward of the center point
        shrink_point = (p + shrink_units * direction).astype(int)
        shrink_points.append(shrink_point)

    s_vals = []
    for p in shrink_points:
        s_vals.append(ros[p[0], p[1]])

    s_vals = np.array(s_vals)
    shrink_points = np.array(shrink_points)

    return shrink_points, s_vals


def simplex_method(step, grow_units=5, shrink_units=5, tol=1e-6, max_iter=1000):
    ## The Simplex method should follow the routine
    ## 1. Generate a random simplex
    ## 2. Compute the area of the inital simplex
    ## 3. Find the smallest value that the simplex holds
    random.seed(10)
    x = np.linspace(-2, 3, step)
    #print(x)
    X, Y = np.meshgrid(x, x)
    ros = rosen([X, Y])
    ## generate 3 random x points (elemnts)
    ## gnerate 3 random y points
    ## Pair them up and create the link
    print(f"ros shape: {ros.shape}")

    random_x = np.random.randint(0 + 200, ros.shape[1] - 500,3)
    random_y = np.random.randint(0 + 200, ros.shape[1] - 500,3)
    points = np.vstack((random_x,random_y)).T

    p_vals = []
    for i in range(points.shape[0]):
        p_vals.append(ros[points[i][0], points[i][1]])
    p_vals = np.array(p_vals)
    paths = []
    paths.append(points)
    for iteration in range(max_iter):
        print(f'Iteration: {iteration}')

        # Sort Points by Function Value
        indices = np.argsort(p_vals)
        points, p_vals = points[indices], p_vals[indices]
        paths.append(points.copy())

        # Perform Reflection
        reflection_points, reflection_vals = reflect(points, p_vals, ros)

        # Evaluate Reflection
        if reflection_vals[0] < p_vals[0]:  # If reflection is better than the best point
            print("Reflection improved.")
            # Perform Expansion (Grow) from the reflection point
            grown_points, g_vals = grow(reflection_points, reflection_vals, grow_units, ros)
            if g_vals[0] < reflection_vals[0]:
                print("Growth is approved.")
                # Use the expanded point if it improves further
                points[-1] = grown_points[0]
                p_vals[-1] = g_vals[0]
            else:
                print("Growth was not approved.")
                # Otherwise, use the reflected point as the updated point
                points[-1] = reflection_points[0]
                p_vals[-1] = reflection_vals[0]
        elif reflection_vals[0] < p_vals[-2]:  # Accept reflection if it's better than the second-worst
            print("Accept reflection if it's better than the second-worst.")
            # Use the reflected point if it shows improvement over the second-worst point
            points[-1] = reflection_points[0]
            p_vals[-1] = reflection_vals[0]
        else:
            print("Performing Shrink.")
            # Perform Shrink: If reflection fails to improve, shrink the simplex toward the best point
            shrink_points, s_vals = shrink(points, p_vals, shrink_units, ros)
            points = shrink_points  # Update all points to the shrunken points
            p_vals = s_vals


        # Check Convergence with Area Tolerance
        area_p = comp_area(points)
        print(f"Area: {area_p}")
        if area_p < tol:
            break

        # Find the minimum point coordinates and value in the function’s domain
    min_indices = points[0]  # This gives us the grid indices of the minimum point
    min_point = (x[min_indices[0]], x[min_indices[1]])  # Convert indices to x, y coordinates
    min_val = p_vals[0]  # Minimum function value at these coordinates

    print(f"Minimum coordinates (x, y): {min_point}")
    print(f"Minimum function value: {min_val}")
    print(f"Iterations: {iteration}")
    return min_point, min_val, paths

if __name__ == "__main__":

    step = 1000
    simplex_method(step)


#
# def plot_animated_simplex_paths(step, n_runs=4, grow_units=5, shrink_units=5, tol=1e-6, max_iter=1000):
#     x = np.linspace(-2, 3, step)
#     X, Y = np.meshgrid(x, x)
#     Z = np.array([[rosen([x_val, y_val]) for x_val in x] for y_val in x])
#
#     fig = plt.figure(figsize=(12, 8))
#     ax = fig.add_subplot(111, projection="3d")
#     ax.plot_surface(X, Y, Z, rstride=5, cstride=5, cmap="viridis", alpha=0.6, edgecolor="none")
#
#     all_paths = []
#     final_points = []
#
#     # Define a fixed set of colors for the simplices
#     simplex_colors = ["blue", "green", "red", "purple", "orange", "brown"]
#
#     def update(i):
#         ax.clear()
#         ax.plot_surface(X, Y, Z, rstride=5, cstride=5, cmap="viridis", alpha=0.6, edgecolor="none")
#
#         # Plot all previous runs with their unique colors
#         for run_index, path in enumerate(all_paths):
#             color = simplex_colors[run_index % len(simplex_colors)]  # Cycle through colors if n_runs > colors list
#             path = np.array(path)
#             for simplex in path:
#                 for j in range(3):
#                     p1, p2 = simplex[j], simplex[(j + 1) % 3]
#                     ax.plot([x[p1[0]], x[p2[0]]], [x[p1[1]], x[p2[1]]],
#                             [rosen([x[p1[0]], x[p1[1]]]), rosen([x[p2[0]], x[p2[1]]])],
#                             color=color, marker="o", markersize=5, alpha=0.6)  # Light color for previous paths
#
#         # Run a new simplex path and plot each step with a unique color
#         if i < n_runs:
#             min_point, min_val, path = simplex_method(step, grow_units, shrink_units, tol, max_iter)
#             final_points.append(min_point)
#             all_paths.append(path)
#             color = simplex_colors[i % len(simplex_colors)]  # Unique color for the current path
#             path = np.array(path)
#             for simplex in path:
#                 for j in range(3):
#                     p1, p2 = simplex[j], simplex[(j + 1) % 3]
#                     ax.plot([x[p1[0]], x[p2[0]]], [x[p1[1]], x[p2[1]]],
#                             [rosen([x[p1[0]], x[p1[1]]]), rosen([x[p2[0]], x[p2[1]]])],
#                             color=color, marker="o", markersize=5, alpha=1.0)  # Color for the active path
#
#             # Plot the final convergence point for the current simplex run
#             ax.scatter(*min_point, min_val, color=color, s=100)
#
#         ax.set_xlabel("X")
#         ax.set_ylabel("Y")
#         ax.set_zlabel("Rosenbrock Value")
#
#     ani = FuncAnimation(fig, update, frames=n_runs, repeat=False)
#     plt.show()
#
# # Run and dynamically plot the animated simplex paths
# step = 1000
# plot_animated_simplex_paths(step)