import csv
import json
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class JSONBoundaryApp:
    
    values_map = {"x": 0, "y": 1, "z": 2}
    reverse_map = {v: k for k, v in values_map.items()}  
    number_of_region = 1
    number_of_elements = 1

    def __init__(self, parent):
        self.selected_index = 0
        self.root = parent
        self._widget_grid_info = {}

        external_time_list = [10, 20, 30, 40, 50]
        external_json_data = {
            "boundary_conditions": {
                "NORTH": {"type": "neumann", "direction":0, "density": 0.0, "reference_position":0.0, "values": [0.0] * 5},
                "SOUTH": {"type": "neumann", "direction":0, "density": 0.0, "reference_position":0.0, "values": [0.0] * 5},
                "WEST": {"type": "neumann", "direction":0, "density": 0.0, "reference_position":0.0, "values": [0.0] * 5},
                "EAST": {"type": "neumann", "direction":0, "density": 0.0, "reference_position":0.0, "values": [0.0] * 5},
                "BOTTOM": {"type": "neumann", "direction":0, "density": 0.0, "reference_position":0.0, "values": [0.0] * 5},
                "TOP": {"type": "neumann", "direction":0, "density": 0.0, "reference_position":0.0, "values": [0.0] * 5},
                "Cavern": {"type": "neumann", "direction":0, "density": 0.0, "reference_position":0.0, "values": [0.0] * 5},
            }
        }

        self.json_data = external_json_data
        self.time_list = external_time_list

        main_frame = tk.Frame(parent)
        main_frame.pack(fill=tk.BOTH, expand=True)

        input_frame = tk.Frame(main_frame)
        input_frame.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.Y)

        plot_frame = tk.Frame(main_frame)
        plot_frame.pack(side=tk.RIGHT, padx=10, pady=10, fill=tk.BOTH, expand=True)

        #self.figure, self.ax = plt.subplots()
        self.figure = Figure()
        self.figure.set_facecolor("#d9d9d9ff")
        self.ax = self.figure.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.enable_drag_and_drop()

        self.dragging_index = None

        # Using grid layout
        self.boundary_label = ttk.Label(input_frame, text="Select Boundary:")
        self.boundary_label.grid(row=0, column=0, sticky="w")
        self.boundary_combobox = ttk.Combobox(input_frame, state="readonly", values=list(self.json_data["boundary_conditions"].keys()))
        self.boundary_combobox.grid(row=0, column=1, sticky="ew")
        self.boundary_combobox.bind("<<ComboboxSelected>>", self.load_boundary_data)

        self.type_label = ttk.Label(input_frame, text="Type:")
        self.type_label.grid(row=1, column=0, sticky="w")
        self.type_combobox = ttk.Combobox(input_frame, values=["neumann", "dirichlet"], state="readonly")
        self.type_combobox.grid(row=1, column=1, sticky="ew")
        self.type_combobox.bind("<<ComboboxSelected>>", self.toggle_fields)

        self.component_label = ttk.Label(input_frame, text="Component/Direction:")
        self.component_label.grid(row=2, column=0, sticky="w")
        self.component_combobox = ttk.Combobox(input_frame, values=["x", "y", "z"], state="readonly")
        self.component_combobox.grid(row=2, column=1, sticky="ew")
        self.component_combobox.bind("<<ComboboxSelected>>", self.update_direction)

        self.density_label = ttk.Label(input_frame, text="Density:")
        self.density_label.grid(row=3, column=0, sticky="w")
        self.density_entry = ttk.Entry(input_frame)
        self.density_entry.grid(row=3, column=1, sticky="ew")
        self.density_entry.bind("<KeyRelease>", self.update_density)

        self.reference_label = ttk.Label(input_frame, text="Ref. Position:")
        self.reference_label.grid(row=4, column=0, sticky="w")
        self.reference_entry = ttk.Entry(input_frame)
        self.reference_entry.grid(row=4, column=1, sticky="ew")
        self.reference_entry.bind("<KeyRelease>", self.update_reference_position)

        self.choose_point_label = ttk.Label(input_frame, text="Choose Point:")
        self.choose_point_label.grid(row=5, column=0, sticky="w")
        self.choose_point_combobox = ttk.Combobox(input_frame, values=self.time_list, state="readonly")
        self.choose_point_combobox.grid(row=5, column=1, sticky="ew")
        self.choose_point_combobox.bind("<<ComboboxSelected>>", self.update_prescribed_value)

        self.pbrowse_button = ttk.Button(input_frame, text="Browse", command=self.browse_csv)
        self.pbrowse_button.grid(row=6, column=1, columnspan=25, pady=5)

        self.prescribed_value_label = ttk.Label(input_frame, text="Prescribed Value:")
        self.prescribed_value_label.grid(row=7, column=0, sticky="w")
        self.prescribed_value_entry = ttk.Entry(input_frame)
        self.prescribed_value_entry.grid(row=7, column=1, sticky="ew")
        self.prescribed_value_entry.bind("<KeyRelease>", self.update_value_and_plot)

        input_frame.columnconfigure(1, weight=1)

        self.boundary_combobox.current(0)
        self.root.bind("<FocusIn>", self.on_focus_in)
        self.load_boundary_data(None)
        for widget in [
            self.density_label, self.density_entry,
            self.reference_label, self.reference_entry
        ]:
            self._widget_grid_info[widget] = widget.grid_info()

    # def SetBoundaryList(self, grid_boundary_names):
    #     #print(type(self.json_data["boundary_conditions"]))
    #     self.json_data["boundary_conditions"] = {key: {"values": [4000000.0] * 5} for key in grid_boundary_names}

    def SetElemNumbers(self,n_region, nu_elem):
        self.number_of_region = n_region
        self.number_of_elements = nu_elem


    def SetBoundaryList(self, grid_boundary_names):
        current_boundaries = self.json_data.get("boundary_conditions", {})
        
        updated_boundaries = {key: value for key, value in current_boundaries.items() if key in grid_boundary_names}
        
        # "NORTH": {"type": "neumann", "direction":0, "density": 0.0, "reference_position":0.0, "values": [4000000.0] * 5},
        for key in grid_boundary_names:
            if key not in updated_boundaries:
                updated_boundaries[key] = {
                        "type": "neumann",
                        "direction": 0,
                        "density": 0.0,
                        "reference_position": 0.0, 
                        "values": [0.0] * len(self.time_list)
                }
        
        self.json_data["boundary_conditions"] = updated_boundaries

        #self.json_data["boundary_conditions"] = grid_boundary_names
        self.update_combobox()
    
    def is_number(self,value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    def is_row_empty(self,row):
        return not row or all(cell.strip() == "" for cell in row)

    def browse_csv(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if file_path:
            try:
                numbers = []
                with open(file_path, newline='', encoding='utf-8') as csvfile:
                    reader = csv.reader(csvfile)
                    for i, row in enumerate(reader):
                        if self.is_row_empty(row):
                            continue
                        if all(self.is_number(cell) for cell in row):
                                numbers.append(float(row[0]))

                    # first_row = next(reader)  # Skip the header row
                    # if all(self.is_number(cell) for cell in first_row):
                    #     try:
                    #         numbers.append(float(first_row[0])) 
                    #     except :
                    #         numbers.append(float(0)) 
                    # for row in reader:
                    #     try:
                    #         numbers.append(float(row[0])) 
                    #     except :
                    #         numbers.append(float(0)) 
                    if len(numbers) != len(self.time_list):
                        messagebox.showerror("Error", f"Time count and Number count are not equal!!!")
                    else:
                        self.update_All_values_and_plot(numbers)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to read file: {str(e)}")

    def update_combobox(self):
        self.boundary_combobox.set("")
        self.boundary_combobox['values'] = list(self.json_data["boundary_conditions"].keys())
        if len(self.boundary_combobox['values']) > 0:
            self.boundary_combobox.current(0)

    def on_focus_in(self, event=None):
        self.enable_drag_and_drop()

    def enable_drag_and_drop(self):
        self.canvas.mpl_connect("motion_notify_event", self.on_drag)
        self.canvas.mpl_connect("button_press_event", self.on_press)
        self.canvas.mpl_connect("button_release_event", self.on_release)

    def update_time_list(self, new_time_list):
        if len(new_time_list) >= len(self.time_list):
            # Add zeros
            for boundary_name, boundary_data in self.json_data["boundary_conditions"].items():
                boundary_data["values"].extend([0.0] * (len(new_time_list) - len(self.time_list)))
        else:
            for boundary_name, boundary_data in self.json_data["boundary_conditions"].items():
                boundary_data["values"] = boundary_data["values"][:(len(new_time_list) - len(self.time_list))] 
        self.time_list = new_time_list
        self.choose_point_combobox["values"] = self.time_list
        if len(self.time_list) > 0:
            self.choose_point_combobox.current(0)
        self.update_prescribed_value()
    
    def update_prescribed_value(self, event=None):
        selected_boundary = self.boundary_combobox.get()
        selected_index = self.choose_point_combobox.current()
        self.selected_index = selected_index  
        if selected_boundary and selected_index >= 0:
            values = self.json_data["boundary_conditions"].get(selected_boundary, {}).get("values", [])
            if selected_index < len(self.time_list):
                if selected_index < len(values):
                    self.prescribed_value_entry.delete(0, tk.END)
                    self.prescribed_value_entry.insert(0, str(values[selected_index]))
                    self.selected_index = selected_index
                else:
                    self.prescribed_value_entry.delete(0, tk.END)
                    self.prescribed_value_entry.insert(0, str("0"))
                    self.json_data["boundary_conditions"][selected_boundary]["values"][selected_index] = 0
                    #self.json_data["boundary_conditions"].get(selected_boundary, {}).get("values", [])
                    self.selected_index = selected_index

            self.update_plot() 

    def load_boundary_data(self, event):
        selected_boundary = self.boundary_combobox.get()
        if selected_boundary:
            boundary_data = self.json_data["boundary_conditions"].get(selected_boundary, {})

            # Set type value
            boundary_type = boundary_data.get("type", "neumann")
            self.type_combobox.set(boundary_type)

            if boundary_type == "neumann":
                # Set direction/component value
                direction = boundary_data.get("direction", 0)  # Default to 1
                self.component_combobox.set(str(self.reverse_map[direction]))

                # Load values
                values = boundary_data.get("values", [])
                ### self.update_time_list(list(range(len(values))))
                self.update_prescribed_value()
                
                # Load other fields (density, reference)
                self.density_entry.delete(0, tk.END)
                self.density_entry.insert(0, boundary_data.get("density", ""))
                
                self.reference_entry.delete(0, tk.END)
                self.reference_entry.insert(0, boundary_data.get("reference_position", ""))
            else:
                # Set direction/component value
                direction = boundary_data.get("component", 0)  # Default to 1
                self.component_combobox.set(str(self.reverse_map[direction]))

                # Load values
                values = boundary_data.get("values", [])
                self.update_prescribed_value()
                
                
        self.update_plot()
        #self.toggle_fields()
        self.toggle_All_widget()

    def update_direction(self, event=None):
        selected_boundary = self.boundary_combobox.get()
        new_direction = self.component_combobox.get()
        selected_type = self.type_combobox.get()
        if selected_boundary and new_direction:
            try:
                if selected_type == "neumann":
                    self.json_data["boundary_conditions"][selected_boundary]["direction"] = int(self.values_map.get(new_direction, 0))
                else:
                    self.json_data["boundary_conditions"][selected_boundary]["component"] = int(self.values_map.get(new_direction, 0))
                self.update_plot()
            except ValueError:
                pass  # Invalid input value

    def update_density(self, event=None):
        selected_boundary = self.boundary_combobox.get()
        new_density = self.density_entry.get()
        if selected_boundary and new_density:
            try:
                self.json_data["boundary_conditions"][selected_boundary]["density"] = float(new_density)
            except:
                self.json_data["boundary_conditions"][selected_boundary]["density"] = float(0)
            self.update_plot()

    def update_reference_position(self, event=None):
        selected_boundary = self.boundary_combobox.get()
        new_reference_position = self.reference_entry.get()
        if selected_boundary and new_reference_position:
            try:
                self.json_data["boundary_conditions"][selected_boundary]["reference_position"] = float(new_reference_position)
            except:
                self.json_data["boundary_conditions"][selected_boundary]["reference_position"] = float(0)
            self.update_plot()

    def update_value_and_plot(self, event=None):
        selected_boundary = self.boundary_combobox.get()
        selected_point = self.choose_point_combobox.current()
        if selected_boundary and selected_point >= 0:
            try:
                new_value = float(self.prescribed_value_entry.get())
                self.json_data["boundary_conditions"][selected_boundary]["values"][selected_point] = new_value
                self.update_plot()
            except ValueError:
                pass  # Invalid input value

    def update_All_values_and_plot(self, numbers):
        selected_boundary = self.boundary_combobox.get()
        if selected_boundary:
            try:
                new_iter = iter(numbers)
                i = 0
                for sp in self.time_list:                
                    self.json_data["boundary_conditions"][selected_boundary]["values"][i] = next(new_iter)
                    i = i + 1
                self.update_plot()
            except ValueError:
                pass  # Invalid input value

    def update_plot(self):
        self.ax.clear()
        selected_boundary = self.boundary_combobox.get()
        if selected_boundary:
            values = self.json_data["boundary_conditions"].get(selected_boundary, {}).get("values", [])
            
            while len(values) < len(self.time_list):
                values.append(0)
            
            self.ax.plot(self.time_list[:len(values)], values[:len(values)], marker='o', linestyle='-', picker=True, label='Values')
            
            if hasattr(self, 'selected_index') and 0 <= self.selected_index < len(self.time_list):
                selected_time = self.time_list[self.selected_index]
                selected_value = values[self.selected_index]
                self.ax.plot(selected_time, selected_value, 'ro', markersize=8, label='Selected Point')

            self.ax.set_xlabel("Time")
            self.ax.set_ylabel("Value")
            self.ax.set_title(f"Boundary: {selected_boundary}")
            self.ax.legend()
            # self.ax.set_facecolor("#d9d9d9ff")
            self.canvas.draw()

    def toggle_All_widget(self):
        selected_type = self.type_combobox.get()
        is_dirichlet = selected_type == "dirichlet"
        def toggle_widget(widget, show):
            if show:
                info = self._widget_grid_info.get(widget)
                if info:
                    widget.grid(**info)
            else:
                widget.grid_remove()

        toggle_widget(self.density_label, not is_dirichlet)
        toggle_widget(self.density_entry, not is_dirichlet)
        toggle_widget(self.reference_label, not is_dirichlet)
        toggle_widget(self.reference_entry, not is_dirichlet)


    def toggle_fields(self, event=None):
        selected_boundary = self.boundary_combobox.get()
        selected_type = self.type_combobox.get()

        if selected_boundary:
            self.json_data["boundary_conditions"][selected_boundary]["type"] = selected_type

        is_dirichlet = selected_type == "dirichlet"
        self.component_label.config(text="Component:" if is_dirichlet else "Direction:")
        if is_dirichlet and selected_boundary:
            bc_data = self.json_data["boundary_conditions"][selected_boundary]
            # remove keys
            bc_data.pop("density", None)
            bc_data.pop("reference_position", None)
            bc_data.pop("direction", None)
            values_tmp = bc_data.pop("values", None)
            
            bc_data["component"] = 0
            bc_data["values"] = values_tmp
        else:
            bc_data = self.json_data["boundary_conditions"][selected_boundary]
            # remove keys
            bc_data.pop( "component", None)
            values_tmp = bc_data.pop("values", None)

            bc_data["direction"] = 0
            bc_data["density"] = 0.0
            bc_data["reference_position"] = 0.0
            bc_data["values"] = values_tmp

        self.toggle_All_widget()


    def on_press(self, event):
        if event.xdata is not None and event.ydata is not None:
            selected_boundary = self.boundary_combobox.get()
            values = self.json_data["boundary_conditions"].get(selected_boundary, {}).get("values", [])

            x_min, x_max = self.ax.get_xlim()
            y_min, y_max = self.ax.get_ylim()

            scale_x = x_max - x_min
            scale_y = y_max - y_min
            
            threshold_x = 0.01 * scale_x  
            threshold_y = 0.01 * scale_y  

            for i, (x, y) in enumerate(zip(self.time_list[:len(values)], values)):
                if abs(event.xdata - x) < threshold_x and abs(event.ydata - y) < threshold_y:
                    self.dragging_index = i
                    break
    def on_drag(self, event):
        if self.dragging_index is not None and event.ydata is not None:
            selected_boundary = self.boundary_combobox.get()
            new_val = float(event.ydata)
            values = self.json_data["boundary_conditions"][selected_boundary]["values"]
            values[self.dragging_index] = new_val

            self.prescribed_value_entry.delete(0, tk.END)
            self.prescribed_value_entry.insert(0, str(new_val))

            self.selected_index = self.dragging_index

            self.update_plot()

    def on_release(self, event):
        self.dragging_index = None


    def on_closing(self):
        plt.close(self.figure)  #   matplotlib
        self.root.quit()
        self.root.destroy()

if __name__ == "__main__":
    root = tk.Tk()
    app = JSONBoundaryApp(root)
    root.protocol("WM_DELETE_WINDOW", app.on_closing)  
    root.mainloop()
