import csv
import json
import os
import sys
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter import scrolledtext
from PIL import Image, ImageTk
import subprocess
import threading
from importlib.resources import files, as_file
from .MyConstitutiveModel import JSONConstitutiveApp
from .MyBoundaryCond import JSONBoundaryApp
from .simulator_runner import SimulatorRunner 
from .script_runner import PythonScriptRunner

# sys.path.append(os.path.join("..", "safeincave"))
from ..Grid import GridHandlerGMSH

sim_in_run = False


def gui():
    # Default file name
    json_file_name = "input_file.json"
    entries = {}
    direction_var = ''
    theta_var = ''
    time_list_text = ''
    grid_boundary_names = {}
    subdomain_region_names = {}
    gn_elems = 0
    time_list_changed = False
    # sim_in_run = False
    list_of_boundary_names = {}
    grid = {}
    number_of_region = 1
    number_of_elements = 1
    current_data = {}
    default_time_list = [10, 20, 30 , 40, 50]

    HERE = os.path.dirname(os.path.abspath(__file__))

    def get_values(entry):
        try:
            text = entry.get().strip()
            if not text:
                return []

            parts = text.split()
            numbers = []

            for part in parts:
                try:
                    num = float(part)
                    if num.is_integer():
                        num = int(num)
                    numbers.append(num)
                except ValueError:
                    pass

            #if len(numbers) == 1:
            #    return numbers[0]
            return numbers
        except Exception as e:
            #print("error in reading input:", e)
            return []
        
    # Populate the form with loaded data
    def populate_form(data):
        # Load Grid settings (path and name)
        if "grid" in data:
            if "name" in data["grid"]:
                grid_name_entry.delete(0, tk.END)
                grid_name_entry.insert(0, str(data["grid"]["name"]))
            if "path" in data["grid"]:
                grid_path_entry.delete(0, tk.END)
                grid_path_entry.insert(0, str(data["grid"]["path"]))

            if "name" in data["grid"] and "path" in data["grid"]:        
                grid_name = "geom"
                #load_grid(grid_name_entry, grid_path_entry, grid_prop_label)

        # Load Output settings (path)
        if "output" in data and "path" in data["output"]:
            output_path_entry.delete(0, tk.END)
            output_path_entry.insert(0, str(data["output"]["path"]))

        # Load simulation settings
        if "simulation_settings" in data:
            for stage, params in data["simulation_settings"].items():
                # Set active state for the stage
                active_var = active_vars[stage]
                active_var.set("Active" if params.get("active", False) else "Inactive")

                # Set hardening for the stage operation
                if stage == "operation":
                    active_var_op = active_vars_operation[stage]
                    active_var_op.set("Active" if params.get("active", False) else "Inactive")

                # Set dt_max with unit conversion
                dt_max = params.get("dt_max", 0)
                time_unit_var = time_unit_vars[stage]
                try:
                    if dt_max < 60:
                        time_unit_var.set("Seconds")
                        entries[f"simulation_settings.{stage}.dt_max"].delete(0, tk.END)
                        entries[f"simulation_settings.{stage}.dt_max"].insert(0, str(dt_max))
                    elif dt_max < 3600:
                        time_unit_var.set("Minutes")
                        entries[f"simulation_settings.{stage}.dt_max"].delete(0, tk.END)
                        entries[f"simulation_settings.{stage}.dt_max"].insert(0, str(dt_max / 60))
                    elif dt_max < 86400:
                        time_unit_var.set("Hours")
                        entries[f"simulation_settings.{stage}.dt_max"].delete(0, tk.END)
                        entries[f"simulation_settings.{stage}.dt_max"].insert(0, str(dt_max / 3600))
                    else:
                        time_unit_var.set("Days")
                        entries[f"simulation_settings.{stage}.dt_max"].delete(0, tk.END)
                        entries[f"simulation_settings.{stage}.dt_max"].insert(0, str(dt_max / 86400))
                except:
                    entries[f"simulation_settings.{stage}.dt_max"].delete(0, tk.END)
                    entries[f"simulation_settings.{stage}.dt_max"].insert(0, str("None"))
                # Set other parameters
                for key, value in params.items():
                    if key not in {"active", "dt_max", "hardening"}:
                        entry_key = f"simulation_settings.{stage}.{key}"
                        entries[entry_key].delete(0, tk.END)
                        entries[entry_key].insert(0, str(value))

        # Load solver settings
        valid_solver_types = {"gKrylovSolver", "LU"}
        if "solver_settings" in data:
            for key, value in data["solver_settings"].items():
                if key == "type":
                    # Set the solver type radio button
                    solver_type_var.set(value if value in valid_solver_types else "KrylovSolver")
                    update_dropdowns()  # Update dependent widgets based on solver type
                elif key == "method":
                    try:
                        # Set the method combobox
                        method_combobox.set(value)
                    except:
                        method_combobox.set("cg")
                elif key == "preconditioner":
                    try:
                        # Set the preconditioner combobox
                        preconditioner_combobox.set(value)
                    except:
                        preconditioner_combobox.set("icc")
                elif key == "relative_tolerance":
                    # Set the relative tolerance entry
                    tolerance_entry.delete(0, tk.END)
                    tolerance_entry.insert(0, str(value))

     
        # Load body force settings
        if "body_force" in data:
            for key, value in data["body_force"].items():
                if key in body_force_parameters:
                    entry_key = f"body_force.{key}"
                    if entry_key in entries:
                        entries[entry_key].delete(0, tk.END)
                        #entries[entry_key].insert(0, str(value))
                        try:
                            entries[entry_key].insert(0, value)
                        except:
                            entries[entry_key].insert(0, str(value))
                elif key == "direction":
                    reverse_mapping = {0: "x", 1: "y", 2: "z"}
                    direction_var.set(reverse_mapping.get(value, "x"))  # Default to "x" if unknown

        # Load time settings
        load_time_settings(data, theta_var, time_list_text)

    # Function to save boundary conditions
    def save_boundary_conditions(data):
        # global boundary_data
        # save_current_boundary_data()  # Save the current boundary's data
        # data["boundary_conditions"] = boundary_data  # Add all boundary data to the JSON structure
        data["boundary_conditions"] = app_bc.json_data["boundary_conditions"]

    # Function to save constitutive model
    def save_constitutive_model(data):
        data["constitutive_model"] = app_cm.data["constitutive_model"]

    def save_to_file(data):
        file_path = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")]
        )
        if file_path:
            with open(file_path, "w") as file:
                json.dump(data, file, indent=4, ensure_ascii=False)
                updateJsonFileName(file_path)
            #messagebox.showinfo("Success", "File saved successfully.")
            display_output("File saved successfully.")

    def scroll_to_line_starting_with(keyword: str):
        lines = left_textbox.get("1.0", tk.END).splitlines()
        for i, line in enumerate(lines, start=1):
            if keyword in line:
                total_lines = int(left_textbox.index("end-1c").split('.')[0])
                position = (i-1) / total_lines
                left_textbox.yview_moveto(position)
                break

    # Function to load an existing JSON file
    def load_from_file():
        global current_data, boundary_data  # Use the global variable
        file_path = filedialog.askopenfilename(
            filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")]
        )
        if file_path:
            if os.path.getsize(file_path) == 0:
                messagebox.showerror("Error", "The selected file is empty!")
                return
            try:
                with open(file_path, "r") as file:
                    current_data = json.load(file)
                    pretty = json.dumps(current_data, indent=4, ensure_ascii=False)
                    left_textbox.delete("1.0", tk.END)
                    left_textbox.insert(tk.END, pretty)

                    populate_form(current_data)

                    if "time_settings" in current_data:
                        if current_data["time_settings"]["time_list"] == "null" or current_data["time_settings"]["time_list"] == None:
                            current_data["time_settings"]["time_list"] = default_time_list
                        # app_bc.time_list = current_data["time_settings"]["time_list"]
                        app_bc.update_time_list(current_data["time_settings"]["time_list"])
                        app_bc.load_boundary_data(None)

                    if "boundary_conditions" in current_data:
                        app_bc.json_data["boundary_conditions"] = current_data["boundary_conditions"]
                        app_bc.update_combobox()
                        app_bc.dragging_index = None
                        app_bc.load_boundary_data(None)

                    if "constitutive_model" in current_data:
                        app_cm.data["constitutive_model"] = current_data["constitutive_model"]
                        app_cm.refresh_tree()

                    updateJsonFileName(file_path)
                
                load_grid(grid_name_entry, grid_path_entry, grid_prop_label)

            except json.JSONDecodeError as e:
                messagebox.showerror("Error", f"Failed to decode the JSON file. Error: {e}")

    # Function to populate time settings
    def load_time_settings(data, theta_var, time_list_text):
        global time_list_changed
        if "time_settings" in data:
            # Load theta
            if "theta" in data["time_settings"]:
                reverse_map = {0:"Fully-implicit", 0.5:"Crank-Nicolson", 1:"Explicit"}
                theta_var.set(reverse_map.get(data["time_settings"]["theta"], "Fully-implicit"))  # Default to "x" if unknown
                #theta_menu.delete(0, tk.END)
                #theta_menu.insert(0, reverse_map.get(data["time_settings"]["theta"], "Fully-implicit"))

                #theta_entry.delete(0, tk.END)
                #theta_entry.insert(0, str(data["time_settings"]["theta"]))

            # Load time list
            if "time_list" in data["time_settings"]:
                times_in_seconds = data["time_settings"]["time_list"]
                times_in_unit = []
                #unit = time_unit_var_TS.get()
                unit = "Seconds"
                try:
                    # Convert times to the selected unit
                    for time in times_in_seconds:
                        if unit == "Seconds":
                            times_in_unit.append(time)
                        elif unit == "Minutes":
                            times_in_unit.append(time / 60)
                        elif unit == "Hours":
                            times_in_unit.append(time / 3600)
                        elif unit == "Days":
                            times_in_unit.append(time / 86400)

                    # Display times in the text box
                    time_list_text.delete(1.0, tk.END)
                    time_list_text.insert(tk.END, "\n".join(map(str, times_in_unit)))
                except:
                    time_list_text.delete(1.0, tk.END)
                time_list_changed = True

    def time_list_key(self, event=None):
        global time_list_changed
        l_data = {}
        l_data["time_settings"] = {}
        l_data["time_settings"]["time_list"] = {}
        save_time_list(l_data)
        time_list_changed = True

    def update_time_list(self, event=None):
        l_data = {}
        l_data["time_settings"] = {}
        l_data["time_settings"]["time_list"] = {}
        save_time_list(l_data)
        app_bc.update_time_list(l_data["time_settings"]["time_list"])
        app_bc.load_boundary_data(None)


    # Function to open file dialog for path selection
    def select_directory(path_entry):
        path = filedialog.askdirectory()
        if path:
            path_entry.delete(0, tk.END)
            path_entry.insert(0, path)

    def load_grid(grid_name_entry, grid_path_entry, grid_prop_label):
        global list_of_boundary_names
        global grid
        global number_of_region
        global number_of_elements
        try:
            # sys.path.append(os.path.join("..", "..", "safeincave"))

            # Get input values
            grid_name = grid_name_entry.get()
            grid_path = grid_path_entry.get()

            # Check if inputs are empty or None
            if not grid_name:
                messagebox.showerror("Input Error", "Grid name must be provided.")
                return
            if not grid_path:
                messagebox.showerror("Input Error", "Grid path must be provided.")
                return

            # Create grid object
            grid = GridHandlerGMSH(grid_name, grid_path)
            list_of_boundary_names = list(grid.get_boundary_names())
            
            number_of_region = grid.n_regions
            number_of_elements = grid.n_elems
            
            app_bc.SetElemNumbers(number_of_region, number_of_elements)
            app_cm.SetElemNumbers(number_of_region, number_of_elements)

            # Extract and display grid information
            grid_boundary_names = grid.get_boundary_names()
            grid_prop_str = "boundary_name: " + str(list(grid_boundary_names)) + "\n"
            app_bc.SetBoundaryList(grid_boundary_names)

            subdomain_region_names = grid.get_subdomain_names()
            grid_prop_str += "subdomain: " + str(list(subdomain_region_names)) + "\n"

            # n_elems = grid.mesh.num_cells()
            # grid_prop_str += "n_elems: " + str(n_elems)
            grid_prop_str += "n_elems: " + str(grid.n_elems)

            grid_prop_label.config(text=grid_prop_str)

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while loading the grid:\n{str(e)}")

    # Save data from the form
    def save_data():
        global grid
        global list_of_boundary_names

        # Check if grid or list_of_boundary_names are not set or empty
        if grid is None or list_of_boundary_names is None or not list_of_boundary_names:
            messagebox.showwarning("Warning", "Grid is not loaded. Attempting to load it now.")
            load_grid(grid_name_entry, grid_path_entry, grid_prop_label)
            return

        data = {}
        try:
            # Save Grid settings
            data["grid"] = {"path": grid_path_entry.get(), "name": grid_name_entry.get()}
            data["output"] = {"path": output_path_entry.get()}

            data["grid"]["regions"] = {value: int(key) for key, value in grid.tags_dict.items()}
            data["grid"]["boundaries"] = list_of_boundary_names

            # Save Solver settings
            if solver_type_var.get() == "KrylovSolver":
                try:
                    tol_f = float(tolerance_entry.get())
                except:
                    tol_f = 0
                data["solver_settings"] = {
                    "type": solver_type_var.get(),
                    "method": method_combobox.get(),
                    "preconditioner": preconditioner_combobox.get(),
                    "relative_tolerance": tol_f
                }
            elif solver_type_var.get() == "LU": #"DirectSolver":
                data["solver_settings"] = {
                    "type": solver_type_var.get(),
                    "method": method_combobox.get(),
                }
            else:
                data["solver_settings"] = {
                    "type": "null",
                    "method": "null",
                    "preconditioner": "null",
                    "relative_tolerance": "null"
                }
                

            # Save Simulation settings
            data["simulation_settings"] = {}
            for stage in {"equilibrium", "operation"}:
                active_var = active_vars[stage]  # Use the dictionary instead of separate variables
                time_unit_var = time_unit_vars[stage]

                # Retrieve dt_max and convert to seconds
                value_str = entries[f"simulation_settings.{stage}.dt_max"].get()

                try:
                    dt_max = float(value_str)
                except:
                    dt_max = None  # نگهداری مقدار خام ورودی (مثل null، رشته خالی یا abc)
                    # dt_max = float(entries[f"simulation_settings.{stage}.dt_max"].get() if entries[f"simulation_settings.{stage}.dt_max"].get() else 0)
                unit = time_unit_var.get()
                if dt_max != None:
                    if unit == "Minutes":
                        dt_max *= 60
                    elif unit == "Hours":
                        dt_max *= 3600
                    elif unit == "Days":
                        dt_max *= 86400

                # Save the stage settings
                data["simulation_settings"][stage] = {
                    "active": active_var.get() == "Active",
                    "dt_max": dt_max,
                }

                # Save additional parameters (e.g., time_tol, n_skip)
                for key in {"time_tol", "n_skip","ite_max"}:
                    entry_key = f"simulation_settings.{stage}.{key}"
                    if entry_key in entries:
                        value = entries[entry_key].get()
                        try:
                            value = eval(value)  # Convert to numeric if possible
                        except:
                            pass
                        data["simulation_settings"][stage][key] = value
                if stage == "operation":
                    key = "hardening"
                    entry_key = f"simulation_settings.{stage}.{key}"
                    value = active_var_op.get() == "Active"
                    data["simulation_settings"][stage][key] = value

            # Save Body Force settings
            data["body_force"] = {}
            for param in body_force_parameters:
                entry_key = f"body_force.{param}"
                if entry_key in entries:
                    if param == "density":
                        value = get_values(entries[entry_key])
                        len_density = len(value)
                        if len_density != 1 and len_density != number_of_region and len_density != number_of_elements:
                            messagebox.showerror("Error", f"Density values length is incorrect!")                
                            #return
                        if len_density == 1:
                            value = value[0]

                    else:
                        value = entries[entry_key].get()
                        try:
                            value = eval(value)  # Convert to numeric if possible
                        except:
                            pass
                       
                    data["body_force"][param] = value

            # Save direction based on dropdown value
            direction_mapping = {"x": 0, "y": 1, "z": 2}
            data["body_force"]["direction"] = int(direction_mapping[direction_var.get()])

            save_time_settings(data)
            save_boundary_conditions(data)
            save_constitutive_model(data)

            save_to_file(data)

            #current_data = json.load(file)
            pretty = json.dumps(data, indent=4, ensure_ascii=False)
            left_textbox.delete("1.0", tk.END)
            left_textbox.insert(tk.END, pretty)
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while saving json:\n{str(e)}")


    # Function to save time settings
    def save_time_settings(data):
        data["time_settings"] = {}

        # Save theta
        values_map = {"Fully-implicit":0, "Crank-Nicolson": 0.5, "Explicit":1}
        data["time_settings"]["theta"] = float(values_map[theta_var.get()] if theta_var.get() else 0)
        
        save_time_list(data)

    def save_time_list(data):
        # Save time list
        raw_times = time_list_text.get(1.0, tk.END).strip().split("\n")
        times_in_seconds = []
        unit = time_unit_var_TS.get()

        # Convert times to seconds
        for raw_time in raw_times:
            try:
                time = float(raw_time)
                if unit == "Minutes":
                    time *= 60
                elif unit == "Hours":
                    time *= 3600
                elif unit == "Days":
                    time *= 86400
                times_in_seconds.append(int(time))
            except ValueError:
                pass  # Ignore invalid inputs

        data["time_settings"]["time_list"] = times_in_seconds

    # Create the main window
    root = tk.Tk()
    root.title("SafeInCave Parameter Manager")
    root.geometry("1900x834")


    # ============ Left Side ============
    left_frame = tk.Frame(root, width=360, height=800, relief=tk.GROOVE, borderwidth=1)
    left_frame.place(x=0, y=34)

    left_textbox = scrolledtext.ScrolledText(left_frame, wrap=tk.WORD)
    left_textbox.place(x=10, y=10, width=340, height=750)

    left_btn1 = ttk.Button(left_frame, text="Load File", command=load_from_file)
    left_btn1.place(x=10, y=765, width=150)

    left_btn2 = ttk.Button(left_frame, text="Save File", command=save_data)
    left_btn2.place(x=180, y=765, width=150)


    # ============ center Frame ============
    center_frame = tk.Frame(root, width=1000, height=800, relief=tk.GROOVE, borderwidth=1)
    center_frame.place(x=360, y=34)
    center_frame.pack_propagate(False)

    # ============ Right Frame ============
    right_frame = tk.Frame(root, width=540, height=800, relief=tk.GROOVE, borderwidth=1)
    right_frame.place(x=1360, y=34)
    right_frame.pack_propagate(False)

    runner = PythonScriptRunner(master=right_frame)
    runner.pack(fill=tk.BOTH, expand=True)

    # Create a notebook (tabs)
    notebook = ttk.Notebook(center_frame)
    notebook.pack(fill=tk.BOTH, expand=True)

    # Tab 1: Grid and Output settings
    grid_output_tab = ttk.Frame(notebook)
    notebook.add(grid_output_tab, text="Grid & Output Settings")

    # Tab 2: Solver Settings
    solver_tab = ttk.Frame(notebook)
    notebook.add(solver_tab, text="Solver Settings")

    # Tab 3: Simulation Settings
    simulation_tab = ttk.Frame(notebook)
    notebook.add(simulation_tab, text="Simulation Settings")

    # Tab 4: Body Force
    body_force_tab = ttk.Frame(notebook)
    notebook.add(body_force_tab, text="Body Force")

    # Tab 5 - Time Settings
    time_settings_tab = ttk.Frame(notebook)
    notebook.add(time_settings_tab, text="Time Settings")

    # # Tab 6 - Boundary Conditions
    # boundary_conditions_tab = ttk.Frame(notebook)
    # notebook.add(boundary_conditions_tab, text="Boundary Conditions")

    # Tab 1 - Grid Name section
    grid_name_frame = ttk.Frame(grid_output_tab)
    grid_name_frame.pack(fill=tk.X, padx=20, pady=10)

    grid_name_label = ttk.Label(grid_name_frame, text="Grid Name:", width=14)
    grid_name_label.pack(side=tk.LEFT, padx=5)

    grid_name_entry = ttk.Entry(grid_name_frame, width=40)
    grid_name_entry.pack(side=tk.LEFT, padx=5)

    # Tab 1 - Grid Path section
    grid_path_frame = ttk.Frame(grid_output_tab)
    grid_path_frame.pack(fill=tk.X, padx=20, pady=10)

    grid_path_label = ttk.Label(grid_path_frame, text="Grid Path:", width=14)
    grid_path_label.pack(side=tk.LEFT, padx=5)

    grid_path_entry = ttk.Entry(grid_path_frame, width=40)
    grid_path_entry.pack(side=tk.LEFT, padx=5)

    grid_path_button = ttk.Button(grid_path_frame, text="Browse", command=lambda: select_directory(grid_path_entry))
    grid_path_button.pack(side=tk.LEFT, padx=5)

    grid_prop_frame = ttk.Frame(grid_output_tab)
    grid_prop_frame.pack(fill=tk.X, padx=20, pady=10)

    grid_prop_label = ttk.Label(grid_prop_frame, text="",wraplength=400, width=56)
    grid_prop_label.pack(side=tk.LEFT, padx=5)

    grid_load_button = ttk.Button(grid_prop_frame, text="Load Grid", command=lambda: load_grid(grid_name_entry, grid_path_entry, grid_prop_label))
    grid_load_button.pack(side=tk.LEFT, padx=5)

    # Tab 1 - Output Path section
    output_frame = ttk.Frame(grid_output_tab)
    output_frame.pack(fill=tk.X, padx=20, pady=10)

    output_path_label = ttk.Label(output_frame, text="Output Path:", width=14)
    output_path_label.pack(side=tk.LEFT, padx=5)

    output_path_entry = ttk.Entry(output_frame, width=40)
    output_path_entry.pack(side=tk.LEFT, padx=5)

    output_path_button = ttk.Button(output_frame, text="Browse", command=lambda: select_directory(output_path_entry))
    output_path_button.pack(side=tk.LEFT, padx=5)


    # Tab 2 - Solver settings section
    solver_frame = ttk.Frame(solver_tab)
    solver_frame.pack(fill=tk.BOTH, padx=20, pady=10)

    solver_parameters = ["type", "method", "preconditioner", "relative_tolerance"]

    solver_type_var = tk.StringVar(value="KrylovSolver")  # Default to KrylovSolver

    method_combobox = None
    preconditioner_combobox = None

    def on_tab_changed(event):
        global time_list_changed
        selected_tab = event.widget.index("current")
        selected_tab_text = notebook.tab(selected_tab, "text")
        json_txt = "grid"
        if selected_tab == 0:
            json_txt = "grid"
        elif selected_tab == 1:
            json_txt = "solver_settings"
        elif selected_tab == 2:
            json_txt = "simulation_settings"
        elif selected_tab == 3:
            json_txt = "body_force"
        elif selected_tab == 4:
            json_txt = "time_settings"
        elif selected_tab == 5:
            json_txt = "boundary_conditions"
        elif selected_tab == 6:
            json_txt = "constitutive_model"
        scroll_to_line_starting_with(json_txt)
        
        if selected_tab_text == "Boundary Conditions" and time_list_changed:
            data = {}
            save_time_settings(data)
            app_bc.update_time_list(data["time_settings"]["time_list"])
            app_bc.load_boundary_data(None)
            time_list_changed = False


    # Function to update dropdown options based on solver type
    def update_dropdowns():
        # Define options based on the solver type
        if solver_type_var.get() == "KrylovSolver":
            method_options = ["cg", "bicg", "bicgstab", "gmres"]
            preconditioner_options = ["icc", "ilu", "petsc_amg", "sor", "hypre", "asm"]
            preconditioner_combobox['state'] = 'readonly'
            tolerance_entry.config(state='normal')
        else:  # DirectSolver
            method_options = ["default", "umfpack", "mumps", "superlu", "petsc"]
            preconditioner_combobox['state'] = 'disabled'
            tolerance_entry.config(state='disabled')
            preconditioner_options = ["icc", "ilu", "petsc_amg", "sor", "hypre"]
        
        # Update method combobox
        method_combobox['values'] = method_options
        method_combobox.set(method_options[0] if method_options else "")
        
        # Update preconditioner combobox
        preconditioner_combobox['values'] = preconditioner_options
        preconditioner_combobox.set(preconditioner_options[0] if preconditioner_options else "")

    # Create Solver Type Section
    type_frame = ttk.Frame(solver_frame)
    type_frame.pack(fill=tk.X, pady=5)

    ttk.Label(type_frame, text="Solver Type:", width=20).pack(side=tk.LEFT)

    ttk.Radiobutton(
        type_frame, text="Krylov Solver", variable=solver_type_var, value="KrylovSolver", command=update_dropdowns
    ).pack(side=tk.LEFT, padx=5)

    ttk.Radiobutton(
        type_frame, text="Direct Solver (LU)", variable=solver_type_var, value="LU", command=update_dropdowns
    ).pack(side=tk.LEFT, padx=5)

    # Create Method Dropdown Section
    method_frame = ttk.Frame(solver_frame)
    method_frame.pack(fill=tk.X, pady=5)

    ttk.Label(method_frame, text="Method:", width=20).pack(side=tk.LEFT)

    method_combobox = ttk.Combobox(method_frame, state="readonly")
    method_combobox.pack(side=tk.LEFT, fill=tk.X, expand=True)

    # Create Preconditioner Dropdown Section
    preconditioner_frame = ttk.Frame(solver_frame)
    preconditioner_frame.pack(fill=tk.X, pady=5)

    ttk.Label(preconditioner_frame, text="Preconditioner:", width=20).pack(side=tk.LEFT)

    preconditioner_combobox = ttk.Combobox(preconditioner_frame, state="readonly")
    preconditioner_combobox.pack(side=tk.LEFT, fill=tk.X, expand=True)

    # Relative tolerance entry
    tolerance_frame = ttk.Frame(solver_frame)
    tolerance_frame.pack(fill=tk.X, pady=5)

    ttk.Label(tolerance_frame, text="Relative Tolerance:", width=20).pack(side=tk.LEFT)
    tolerance_entry = ttk.Entry(tolerance_frame)
    tolerance_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)

    # Initialize the dropdowns
    update_dropdowns()

    # Tab 3 - Simulation settings section
    simulation_frame = ttk.Frame(simulation_tab)
    simulation_frame.pack(fill=tk.BOTH, expand=True)

    canvas = tk.Canvas(simulation_frame)
    scrollbar = ttk.Scrollbar(simulation_frame, orient=tk.VERTICAL, command=canvas.yview)
    scrollable_frame = ttk.Frame(canvas)

    scrollable_frame.bind(
        "<Configure>",
        lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
    )

    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)

    canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    # Shared active variables for equilibrium and operation
    active_vars = {
        "equilibrium": tk.StringVar(value="Active"),
        "operation": tk.StringVar(value="Active"),
    }
    active_vars_operation = {
        "operation": tk.StringVar(value="Active"),
    }

    time_unit_vars = {
        "equilibrium": tk.StringVar(value="Seconds"),
        "operation": tk.StringVar(value="Seconds"),
    }

    # Define simulation settings structure
    simulation_parameters = {
        "equilibrium": ["time_tol"],
        "operation": ["n_skip"]
    }

    # Generate simulation settings dynamically
    for stage, additional_params in simulation_parameters.items():
        section_label = ttk.Label(scrollable_frame, text=f"{stage.capitalize()} Stage", font=("Arial", 12, "bold"))
        section_label.pack(anchor="w", pady=10)

        # Active/Inactive radio buttons
        active_frame = ttk.Frame(scrollable_frame)
        active_frame.pack(fill=tk.X, padx=20, pady=5)

        active_var = active_vars[stage]
        ttk.Label(active_frame, text="State:", width=20).pack(side=tk.LEFT)

        ttk.Radiobutton(active_frame, text="Active", variable=active_var, value="Active").pack(side=tk.LEFT)
        ttk.Radiobutton(active_frame, text="Inactive", variable=active_var, value="Inactive").pack(side=tk.LEFT)

        # Time Step Size (dt_max) with unit selection
        time_step_frame = ttk.Frame(scrollable_frame)
        time_step_frame.pack(fill=tk.X, padx=20, pady=5)

        ttk.Label(time_step_frame, text="Time Step Size:", width=20).pack(side=tk.LEFT)

        dt_max_entry = ttk.Entry(time_step_frame)
        dt_max_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        entries[f"simulation_settings.{stage}.dt_max"] = dt_max_entry

        time_unit_var = time_unit_vars[stage]
        time_unit_menu = ttk.OptionMenu(time_step_frame, time_unit_var, "Seconds", "Seconds", "Minutes", "Hours", "Days")
        time_unit_menu.pack(side=tk.LEFT)

        if stage == "equilibrium":
            param = "ite_max"
            param_frame = ttk.Frame(scrollable_frame)
            param_frame.pack(fill=tk.X, padx=20, pady=5)

            ttk.Label(param_frame, text=param, width=20).pack(side=tk.LEFT)

            param_entry = ttk.Entry(param_frame)
            param_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
            entries[f"simulation_settings.{stage}.{param}"] = param_entry

        # Additional parameters
        for param in additional_params:
            param_frame = ttk.Frame(scrollable_frame)
            param_frame.pack(fill=tk.X, padx=20, pady=5)

            ttk.Label(param_frame, text=param, width=20).pack(side=tk.LEFT)

            param_entry = ttk.Entry(param_frame)
            param_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
            entries[f"simulation_settings.{stage}.{param}"] = param_entry
        if stage == "operation":
            # Active/Inactive radio buttons
            param = "hardening"
            active_frame = ttk.Frame(scrollable_frame)
            active_frame.pack(fill=tk.X, padx=20, pady=5)

            active_var_op = active_vars_operation[stage]
            ttk.Label(active_frame, text="hardening:", width=20).pack(side=tk.LEFT)

            ttk.Radiobutton(active_frame, text="Active", variable=active_var_op, value="Active").pack(side=tk.LEFT)
            ttk.Radiobutton(active_frame, text="Inactive", variable=active_var_op, value="Inactive").pack(side=tk.LEFT)
            #entries[f"simulation_settings.{stage}.{param}"] = param_entry
            


    # Tab 4 - Body Force section
    body_force_frame = ttk.Frame(body_force_tab)
    body_force_frame.pack(fill=tk.BOTH, padx=20, pady=10)

    # Body force parameters
    body_force_parameters = ["gravity", "density"]

    # Create entry fields for gravity and density
    for param in body_force_parameters:
        param_frame = ttk.Frame(body_force_frame)
        param_frame.pack(fill=tk.X, pady=5)

        if param == "gravity":
            ttk.Label(param_frame, text="Gravity (m/s²):", width=20).pack(side=tk.LEFT)
        if param == "density":
            ttk.Label(param_frame, text="Density (kg/m³):", width=20).pack(side=tk.LEFT)

        param_entry = ttk.Entry(param_frame)
        param_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        entries[f"body_force.{param}"] = param_entry

    # Dropdown for direction
    direction_frame = ttk.Frame(body_force_frame)
    direction_frame.pack(fill=tk.X, pady=5)

    ttk.Label(direction_frame, text="Direction:", width=20).pack(side=tk.LEFT)

    direction_var = tk.StringVar(value="x")  # Default value is "x"
    direction_menu = ttk.OptionMenu(direction_frame, direction_var, "x", "x", "y", "z")
    direction_menu.pack(side=tk.LEFT, fill=tk.X, expand=True)

    # Tab 5 - Time Settings section
    time_settings_frame = ttk.Frame(time_settings_tab)
    time_settings_frame.pack(fill=tk.BOTH, padx=20, pady=10)

    # Theta parameter
    theta_frame = ttk.Frame(time_settings_frame)
    theta_frame.pack(fill=tk.X, pady=5)

    ttk.Label(theta_frame, text="Time integration scheme: ", width=30).pack(side=tk.LEFT)

    theta_var = tk.StringVar(value="x")  # Default value is "x"
    theta_menu = ttk.OptionMenu(theta_frame, theta_var, "Fully-implicit", "Fully-implicit", "Crank-Nicolson", "Explicit")
    theta_menu.pack(side=tk.LEFT, fill=tk.X, expand=True)

    # Time List parameter
    time_list_frame = ttk.Frame(time_settings_frame)
    time_list_frame.pack(fill=tk.BOTH, pady=5)

    ttk.Label(time_list_frame, text="Time List:", width=20).pack(side=tk.LEFT, anchor="n")

    time_list_text = tk.Text(time_list_frame, height=5, width=40)
    time_list_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)
    time_list_text.bind("<KeyRelease>", time_list_key)

    time_list_text.delete(1.0, tk.END)
    time_list_text.insert(tk.END, "\n".join(map(str, [10, 20, 30 , 40, 50])))
    current_data["time_settings"] = {}
    current_data["time_settings"]["time_list"] = [10, 20, 30 , 40, 50]

    time_unit_var_TS = tk.StringVar(value="Seconds")  # Default unit is Seconds
    time_unit_menu = ttk.OptionMenu(time_list_frame, time_unit_var_TS, "Seconds", "Seconds", "Minutes", "Hours", "Days")
    time_unit_menu.pack(side=tk.LEFT, padx=5, anchor="n")

    timebrowse_button = ttk.Button(time_list_frame, text="Browse", command=lambda: browse_csv_time())
    timebrowse_button.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5, anchor="s")


    Boundary_Cond_tab = ttk.Frame(notebook)
    notebook.add(Boundary_Cond_tab, text="Boundary Conditions")

    Constitutive_model_tab = ttk.Frame(notebook)
    notebook.add(Constitutive_model_tab, text="Constitutive model")

    notebook.bind("<<NotebookTabChanged>>", on_tab_changed)  


    app_cm = JSONConstitutiveApp(Constitutive_model_tab)
    app_bc = JSONBoundaryApp(Boundary_Cond_tab)



    # logo_image = Image.open(os.path.join(HERE, "logo_2.png")).convert("RGBA") 
    with as_file(files("safeincave.app.assets") / "logo_2.png") as logo_2_path:
        logo_image = Image.open(logo_2_path).convert("RGBA")
        logo_image = logo_image.resize((190, 32)) 
        logo_photo = ImageTk.PhotoImage(logo_image)

    logo_label = tk.Label(root, image=logo_photo)
    # logo_label.place(x=960-105, y=0)
    logo_label.place(x=0, y=0)

    #root.iconphoto(False, logo_photo)
    with as_file(files("safeincave.app.assets") / "logo_alone_2.png") as logo32_path:
        # logo_image = Image.open(logo32_path).convert("RGBA")
        # icon = tk.PhotoImage(logo_image)
        icon = tk.PhotoImage(file=str(logo32_path))
        root.iconphoto(False, icon)


    style = ttk.Style()
    style.configure("TNotebook.Tab", padding=[5, 8])

    # Buttons for file operations
    button_frame = ttk.Frame(center_frame)
    button_frame.pack(fill=tk.X, padx=10, pady=10)

    # load_button = ttk.Button(button_frame, text="Load File", command=load_from_file)
    # load_button.pack(side=tk.LEFT, padx=5)

    # save_button = ttk.Button(button_frame, text="Save File", command=save_data)
    # save_button.pack(side=tk.LEFT, padx=5)

    text_widget = scrolledtext.ScrolledText(center_frame, wrap=tk.WORD)
    text_widget.pack(padx=5, pady=5, fill="both", expand=True)  

    def display_output(text):
        text_widget.insert(tk.END, text)
        text_widget.see(tk.END) 

    def run_simulation():
        global sim_in_run
        if sim_in_run:
            sim_runner.stop()
            sim_in_run = False
            btn_run.config(text="Run Simulation")
        else:
            file_path = filedialog.askopenfilename(
                filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")]
            )
            if file_path:
                updateJsonFileName(file_path)
                sim_runner.setJsonFile(file_path)
                btn_run.config(text="Stop Simulation")
                sim_in_run = True
                sim_runner.run()


    def browse_csv_time():
        global time_list_changed
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if file_path:
            try:
                time_numbers = []
                with open(file_path, newline='', encoding='utf-8') as csvfile:
                    reader = csv.reader(csvfile)
                    next(reader)  # Skip the header row
                    for row in reader:
                        time_numbers.append(row[0])  
                    
                    time_list_text.delete(1.0, tk.END)
                    time_list_text.insert(tk.END, "\n".join(map(str, time_numbers)))
                    time_list_changed = True
                    update_time_list(None)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to read file: {str(e)}")

    sim_runner = SimulatorRunner(display_output)


    ttk.Label(button_frame, text="JSON File:", width=10).pack(side=tk.LEFT)
    label_file = ttk.Label(button_frame, text=json_file_name, width=50)
    label_file.pack(side=tk.LEFT)


            
    def updateJsonFileName(filename1):
        json_file_name = os.path.basename(filename1)
        #json_file_name = filename1
        label_file.config(text=json_file_name)

    #exit_button = ttk.Button(button_frame, text="Exit", command=root.destroy)
    #exit_button.pack(side=tk.RIGHT, padx=5)

    btn_run = tk.Button(button_frame, text="Run Simulation", command=run_simulation)
    btn_run.pack(side=tk.RIGHT, padx=5)

    # Start the GUI loop

    root.protocol("WM_DELETE_WINDOW", app_bc.on_closing) 
    root.mainloop()

if __name__ == '__main__':
    gui()