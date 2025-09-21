import json
import csv
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import locale

import locale
locale.setlocale(locale.LC_ALL, '')  # یعنی استفاده از locale پیش‌فرض سیستم


class JSONConstitutiveApp:
    in_edit_mode = False
    number_of_region = 1
    number_of_elements = 1
    
    def __init__(self, parent):
        self.root = parent
        #self.root.title("Constitutive model")
        self.data = {"constitutive_model": {"elastic": {}, "nonelastic": {}}}
        self.entries = []        
        self.vars = []
        self.btn_style = ttk.Style()
        self.btn_style.configure("Red.TButton", background="red")
        self.loading = True  # فلگ برای تشخیص لود اولیه

        self.create_widgets()

    def SetElemNumbers(self,n_region, nu_elem):
        self.number_of_region = n_region
        self.number_of_elements = nu_elem

    def is_number(self,value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    def is_row_empty(self,row):
        return not row or all(cell.strip() == "" for cell in row)

    def browse_csv(self,entry_widget):
        self.loading = True
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
                                numbers.append(row[0]) 
                    entry_widget.delete(0, tk.END)
                    entry_widget.insert(0, " ".join(numbers))
            except Exception as e:
                messagebox.showerror("Error", f"Failed to read file: {str(e)}")
        self.loading = False

    def on_var_change(self, entry):
        if self.loading:
            return        
        entry.config(bg='orange')
        self.add_button.configure(style="Red.TButton")

    def reset_colors(self):
        for entry in self.entries:
            entry.config(bg='white')
        self.add_button.configure(style="TButton") 

    def _make_entry(self):
        var = tk.StringVar()
        entry = tk.Entry(self.root, textvariable=var)
        var.trace_add('write', lambda var_name, index, mode, e=entry: self.on_var_change(e))
        self.entries.append(entry)
        self.vars.append(var)
        return entry
    
    def create_widgets(self):

        self.type_label = ttk.Label(self.root, text="Select Type:")
        self.type_label.grid(row=0, column=0, padx=5, pady=5, sticky="w")

        self.type_combo = ttk.Combobox(self.root, values=["Spring", "KelvinVoigt", "ViscoplasticDesai", "DislocationCreep"], state="readonly")
        self.type_combo.grid(row=0, column=1, padx=5, pady=5, sticky="w")
        self.type_combo.bind("<<ComboboxSelected>>", self.on_type_selected)
        
        self.add_button = ttk.Button(self.root, text="Add Model", command=self.add_model)
        self.add_button.grid(row=1, column=3, padx=5, pady=5, sticky="w")


        self.equi_var = tk.BooleanVar(value=False)
        self.equ_label = tk.Label(self.root, text="equilibrium :")
        self.radio_true = tk.Radiobutton(self.root, text="True", variable=self.equi_var, value=True)
        self.radio_false = tk.Radiobutton(self.root, text="False", variable=self.equi_var, value=False)



        self.e_entry = self._make_entry()
        self.nu_entry = self._make_entry()
        self.eta_entry = self._make_entry()
        self.a_entry = self._make_entry()
        self.n_entry = self._make_entry()
        self.q_entry = self._make_entry()
        self.t_entry = self._make_entry()
        self.mu_1_entry = self._make_entry()
        self.N_1_entry = self._make_entry()
        self.n_desai_entry = self._make_entry()
        self.a_1_entry = self._make_entry()
        self.eta_desai_entry = self._make_entry()
        self.beta_1_entry = self._make_entry()
        self.beta_entry = self._make_entry()
        self.m_entry = self._make_entry()
        self.gamma_entry = self._make_entry()
        self.alpha_0_entry = self._make_entry()
        self.k_v_entry = self._make_entry()
        self.sigma_t_entry = self._make_entry()


        self.e_label = ttk.Label(self.root, text="E:")
        #self.e_entry = ttk.Entry(self.root)
        self.e_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.e_entry))

        self.nu_label = ttk.Label(self.root, text="nu:")
        #self.nu_entry = ttk.Entry(self.root)
        self.nu_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.nu_entry))

        self.eta_label = ttk.Label(self.root, text="eta:")
        #self.eta_entry = ttk.Entry(self.root)
        self.eta_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.eta_entry))

        self.a_label = ttk.Label(self.root, text="A:")
        #self.a_entry = ttk.Entry(self.root)
        self.a_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.a_entry))

        self.n_label = ttk.Label(self.root, text="n:")
        #self.n_entry = ttk.Entry(self.root)
        self.n_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.n_entry))

        self.q_label = ttk.Label(self.root, text="Q:")
        #self.q_entry = ttk.Entry(self.root)
        self.q_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.q_entry))

        self.t_label = ttk.Label(self.root, text="T:")
        #self.t_entry = ttk.Entry(self.root)
        self.t_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.t_entry))

        # Labels and Entry fields for variables
        self.mu_1_label = ttk.Label(self.root, text="μ₁:")
        #self.mu_1_entry = ttk.Entry(self.root)
        self.mu_1_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.mu_1_entry))

        self.N_1_label = ttk.Label(self.root, text="N₁:")
        #self.N_1_entry = ttk.Entry(self.root)
        self.N_1_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.N_1_entry))

        self.n_desai_label = ttk.Label(self.root, text="n:")
        #self.n_desai_entry = ttk.Entry(self.root)
        self.n_desai_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.n_desai_entry))

        self.a_1_label = ttk.Label(self.root, text="a₁:")
        #self.a_1_entry = ttk.Entry(self.root)
        self.a_1_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.a_1_entry))

        self.eta_desai_label = ttk.Label(self.root, text="η:")
        #self.eta_desai_entry = ttk.Entry(self.root)
        self.eta_desai_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.eta_desai_entry))

        self.beta_1_label = ttk.Label(self.root, text="β₁:")
        #self.beta_1_entry = ttk.Entry(self.root)
        self.beta_1_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.beta_1_entry))

        self.beta_label = ttk.Label(self.root, text="β:")
        #self.beta_entry = ttk.Entry(self.root)
        self.beta_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.beta_entry))

        self.m_label = ttk.Label(self.root, text="m:")
        #self.m_entry = ttk.Entry(self.root)
        self.m_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.m_entry))

        self.gamma_label = ttk.Label(self.root, text="γ:")
        #self.gamma_entry = ttk.Entry(self.root)
        self.gamma_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.gamma_entry))

        self.alpha_0_label = ttk.Label(self.root, text="α₀:")
        #self.alpha_0_entry = ttk.Entry(self.root)
        self.alpha_0_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.alpha_0_entry))

        self.k_v_label = ttk.Label(self.root, text="k_v:")
        #self.k_v_entry = ttk.Entry(self.root)
        self.k_v_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.k_v_entry))

        self.sigma_t_label = ttk.Label(self.root, text="σ_t:")
        #self.sigma_t_entry = ttk.Entry(self.root)
        self.sigma_t_browse_button = tk.Button(self.root, text="Browse", command=lambda: self.browse_csv(self.sigma_t_entry))



        self.name_label = ttk.Label(self.root, text="Choose Name:")
        self.name_label.grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.name_entry = ttk.Entry(self.root)
        self.name_entry.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        self.tree = ttk.Treeview(self.root, columns=("Category", "Name", "Type", "Edit", "Remove"), show="headings")
        self.tree.heading("Category", text="Category")
        self.tree.heading("Name", text="Name")
        self.tree.heading("Type", text="Type")
        self.tree.heading("Edit", text="Edit")
        self.tree.heading("Remove", text="Remove")
        self.tree.grid(row=2, column=0, rowspan=6, columnspan=6, sticky="nsew", padx=5, pady=5)
        self.tree.bind('<Button-1>', self.handle_action)
        self.tree.column("Category", width=80, minwidth=50, stretch=False)
        self.tree.column("Name", width=100, minwidth=50, stretch=False)
        self.tree.column("Type", width=100, minwidth=50, stretch=False)
        self.tree.column("Edit", width=60, minwidth=50, stretch=False)
        self.tree.column("Remove", width=60, minwidth=50, stretch=False)

        self.add_button = ttk.Button(self.root, text="Apply Changes", command=self.keep_changes)
        self.add_button.grid(row=10, column=2, padx=5, pady=5, sticky="w")



        # self.save_button = ttk.Button(self.root, text="Save to File", command=self.save_to_file)
        # self.save_button.grid(row=10, column=0, padx=5, pady=5, sticky="w")
        # self.load_button = ttk.Button(self.root, text="Load from File", command=self.load_from_file)
        # self.load_button.grid(row=10, column=1, padx=5, pady=5, sticky="w")

        self.refresh_tree()

    def on_type_selected(self, event):
        self.type_select_change()

    def add_entry_with_button(self, label, entry, button, row, column_label, column_entry, column_button):
        label.grid(row=row, column=column_label, padx=5, pady=5, sticky="w")
        entry.grid(row=row, column=column_entry, padx=5, pady=5, sticky="w")
        if button:  # Check if button is not None
            button.grid(row=row, column=column_button, padx=5, pady=5, sticky="w")

    def type_select_change(self):
        selected_type = self.type_combo.get()
        self.hide_all_parameters()

        self.equ_label.grid(row=0, column=7, padx=5, pady=5, sticky='w')
        self.radio_true.grid(row=0, column=8, padx=5, pady=5)
        self.radio_false.grid(row=0, column=9, padx=5, pady=5)


        if selected_type == "Spring":
            self.add_entry_with_button(self.e_label, self.e_entry, self.e_browse_button, 1, 7, 8, 9)
            self.add_entry_with_button(self.nu_label, self.nu_entry, self.nu_browse_button, 2, 7, 8, 9)
        elif selected_type == "KelvinVoigt":
            self.add_entry_with_button(self.e_label, self.e_entry, self.e_browse_button, 1, 7, 8, 9)
            self.add_entry_with_button(self.nu_label, self.nu_entry, self.nu_browse_button, 2, 7, 8, 9)
            self.add_entry_with_button(self.eta_label, self.eta_entry, self.eta_browse_button, 3, 7, 8, 9)
        elif selected_type == "DislocationCreep":
            self.add_entry_with_button(self.a_label, self.a_entry, self.a_browse_button, 1, 7, 8, 9)
            self.add_entry_with_button(self.n_label, self.n_entry, self.n_browse_button, 2, 7, 8, 9)
            self.add_entry_with_button(self.q_label, self.q_entry, self.q_browse_button, 3, 7, 8, 9)
            self.add_entry_with_button(self.t_label, self.t_entry, self.t_browse_button, 4, 7, 8, 9)
        elif selected_type == "ViscoplasticDesai":
            self.add_entry_with_button(self.mu_1_label, self.mu_1_entry, self.mu_1_browse_button, 1, 7, 8, 9)
            self.add_entry_with_button(self.N_1_label, self.N_1_entry, self.N_1_browse_button, 2, 7, 8, 9)
            self.add_entry_with_button(self.n_desai_label, self.n_desai_entry, self.n_desai_browse_button, 3, 7, 8, 9)
            self.add_entry_with_button(self.a_1_label, self.a_1_entry, self.a_1_browse_button, 4, 7, 8, 9)
            self.add_entry_with_button(self.eta_desai_label, self.eta_desai_entry, self.eta_desai_browse_button, 5, 7, 8, 9)
            self.add_entry_with_button(self.beta_1_label, self.beta_1_entry, self.beta_1_browse_button, 6, 7, 8, 9)
            self.add_entry_with_button(self.beta_label, self.beta_entry, self.beta_browse_button, 7, 7, 8, 9)
            self.add_entry_with_button(self.m_label, self.m_entry, self.m_browse_button, 8, 7, 8, 9)
            self.add_entry_with_button(self.gamma_label, self.gamma_entry, self.gamma_browse_button, 9, 7, 8, 9)
            self.add_entry_with_button(self.alpha_0_label, self.alpha_0_entry, self.alpha_0_browse_button, 10, 7, 8, 9)
            self.add_entry_with_button(self.k_v_label, self.k_v_entry, self.k_v_browse_button, 11, 7, 8, 9)
            self.add_entry_with_button(self.sigma_t_label, self.sigma_t_entry, self.sigma_t_browse_button, 12, 7, 8, 9)
    
    def hide_all_parameters(self):
        self.e_label.grid_remove()
        self.e_entry.grid_remove()
        self.e_browse_button.grid_remove()
        self.nu_label.grid_remove()
        self.nu_entry.grid_remove()
        self.eta_label.grid_remove()
        self.eta_entry.grid_remove()
        self.a_label.grid_remove()
        self.a_entry.grid_remove()
        self.n_label.grid_remove()
        self.n_entry.grid_remove()
        self.q_label.grid_remove()
        self.q_entry.grid_remove()
        self.t_label.grid_remove()
        self.t_entry.grid_remove()
        # Hiding Labels and Entry fields using grid_remove
        self.mu_1_label.grid_remove()
        self.mu_1_entry.grid_remove()
        self.N_1_label.grid_remove()
        self.N_1_entry.grid_remove()
        self.n_desai_label.grid_remove()
        self.n_desai_entry.grid_remove()
        self.a_1_label.grid_remove()
        self.a_1_entry.grid_remove()
        self.eta_desai_label.grid_remove()
        self.eta_desai_entry.grid_remove()
        self.beta_1_label.grid_remove()
        self.beta_1_entry.grid_remove()
        self.beta_label.grid_remove()
        self.beta_entry.grid_remove()
        self.m_label.grid_remove()
        self.m_entry.grid_remove()
        self.gamma_label.grid_remove()
        self.gamma_entry.grid_remove()
        self.alpha_0_label.grid_remove()
        self.alpha_0_entry.grid_remove()
        self.k_v_label.grid_remove()
        self.k_v_entry.grid_remove()
        self.sigma_t_label.grid_remove()
        self.sigma_t_entry.grid_remove()

        self.e_browse_button.grid_remove()
        self.nu_browse_button.grid_remove()
        self.eta_browse_button.grid_remove()
        self.a_browse_button.grid_remove()
        self.n_browse_button.grid_remove()
        self.q_browse_button.grid_remove()
        self.t_browse_button.grid_remove()
        self.mu_1_browse_button.grid_remove()
        self.N_1_browse_button.grid_remove()
        self.n_desai_browse_button.grid_remove()
        self.a_1_browse_button.grid_remove()
        self.eta_desai_browse_button.grid_remove()
        self.beta_1_browse_button.grid_remove()
        self.beta_browse_button.grid_remove()
        self.m_browse_button.grid_remove()
        self.gamma_browse_button.grid_remove()
        self.alpha_0_browse_button.grid_remove()
        self.k_v_browse_button.grid_remove()
        self.sigma_t_browse_button.grid_remove()

    def add_model(self):
        selected_type = self.type_combo.get()
        model_name = self.name_entry.get()
        
        return self.Add_Keep_Changes(selected_type, model_name, model_name)

    def Add_Keep_Changes(self, selected_type, model_name, new_model_name):
        if model_name and selected_type:
            parameters = {}
            
            def get_values(entry,param_name):
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

                    len_ = len(numbers)
                    if len_ != 1 and len_ != self.number_of_region and len_ != self.number_of_elements:
                        messagebox.showerror("Warning", f"{param_name} values length is incorrect!")                
                    if len(numbers) == 1:
                        return numbers[0]
                    return numbers
                except Exception as e:
                    #print("error in reading input:", e)
                    return []

            #model_name = 
            category1 = "elastic"
            if selected_type in ["Spring", "KelvinVoigt"]:
                category1 = "elastic"
                parameters["E"] = get_values(self.e_entry,"E")
                parameters["nu"] = get_values(self.nu_entry,"nu")
            
            if selected_type == "KelvinVoigt":
                category1 = "nonelastic"
                parameters["eta"] = get_values(self.eta_entry,"eta")
            
            if selected_type == "DislocationCreep":
                category1 = "nonelastic"
                parameters["A"] = get_values(self.a_entry,"A")
                parameters["n"] = get_values(self.n_entry,"n")
                parameters["Q"] = get_values(self.q_entry,"Q")
                parameters["T"] = get_values(self.t_entry,"T")
            
            if selected_type == "ViscoplasticDesai":
                category1 = "nonelastic"
                parameters["mu_1"]    = get_values(self.mu_1_entry,    "mu_1")
                parameters["N_1"]     = get_values(self.N_1_entry,     "N_1")
                parameters["n"]       = get_values(self.n_desai_entry, "n")
                parameters["a_1"]     = get_values(self.a_1_entry,     "a_1")
                parameters["eta"]     = get_values(self.eta_desai_entry, "eta")
                parameters["beta_1"]  = get_values(self.beta_1_entry,  "beta_1")
                parameters["beta"]    = get_values(self.beta_entry,    "beta")
                parameters["m"]       = get_values(self.m_entry,       "m")
                parameters["gamma"]   = get_values(self.gamma_entry,   "gamma")
                parameters["alpha_0"] = get_values(self.alpha_0_entry, "alpha_0")
                parameters["k_v"]     = get_values(self.k_v_entry,     "k_v")
                parameters["sigma_t"] = get_values(self.sigma_t_entry, "sigma_t")

            if model_name != new_model_name:
                self.data["constitutive_model"][category1][new_model_name] = self.data["constitutive_model"][category1].pop(model_name)

            self.data["constitutive_model"][category1][new_model_name] = {
                "type": selected_type,
                "active": True,
                "equilibrium": self.equi_var.get(),
                "parameters": parameters
            }
            
            self.refresh_tree()
        else:
            messagebox.showwarning("Input Error", "Please provide a name and select a type.")
        self.reset_colors()

    def refresh_tree(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        for category, models in self.data["constitutive_model"].items():
            for name, properties in models.items():
                self.tree.insert("", tk.END, values=(category, name, properties.get("type", ""), "Edit", "Remove"))

    item_keep_name = ""
    item_keep_type = ""
    def handle_action(self, event):
        region = self.tree.identify_region(event.x, event.y)
        column = self.tree.identify_column(event.x)
        item_id = self.tree.identify_row(event.y)
        self.in_edit_mode = False
        if region == 'cell' and item_id:
            values = self.tree.item(item_id, 'values')
            if column == "#4":
                self.edit_item_popup(values[0], values[1])
                self.in_edit_mode = True
                self.item_keep_name = values[1]
                self.item_keep_type = values[2]

            elif column == "#5":
                self.confirm_and_remove(values[0], values[1])
            else:
                self.hide_all_parameters()
    
    def keep_changes(self):
        if self.in_edit_mode:
            if(self.confirm_keep_item_changes(self.item_keep_name,self.item_keep_type)):
                new_model_name = self.name_entry.get()
                self.Add_Keep_Changes(self.item_keep_type, self.item_keep_name, new_model_name)
            
    def edit_item_popup(self, category, name):
        self.reset_colors()
        self.loading = True
        model = self.data["constitutive_model"][category].get(name, {})
        self.type_combo.set(model.get("type", ""))
        self.name_entry.delete(0, tk.END)
        self.name_entry.insert(0, name)
        self.type_select_change()
        self.equi_var.set(model.get("equilibrium", False))
        # Show parameters only for Spring type
        if model.get("type") == "Spring":
            self.e_entry.delete(0, tk.END)
            self.e_entry.insert(0, model["parameters"].get("E", []) if "E" in model["parameters"] else "")
            self.nu_entry.delete(0, tk.END)
            self.nu_entry.insert(0, model["parameters"].get("nu", []) if "nu" in model["parameters"] else "")
        elif model.get("type") == "KelvinVoigt":
            self.e_entry.delete(0, tk.END)
            self.e_entry.insert(0, model["parameters"].get("E", []) if "E" in model["parameters"] else "")
            self.nu_entry.delete(0, tk.END)
            self.nu_entry.insert(0, model["parameters"].get("nu", []) if "nu" in model["parameters"] else "")
            self.eta_entry.delete(0, tk.END)
            self.eta_entry.insert(0, model["parameters"].get("eta", []) if "nu" in model["parameters"] else "")
        elif model.get("type") == "DislocationCreep":
            self.a_entry.delete(0, tk.END)
            self.a_entry.insert(0, model["parameters"].get("A", []) if "A" in model["parameters"] else "")

            self.n_entry.delete(0, tk.END)
            self.n_entry.insert(0, model["parameters"].get("n", []) if "n" in model["parameters"] else "")

            self.q_entry.delete(0, tk.END)
            self.q_entry.insert(0, model["parameters"].get("Q", []) if "Q" in model["parameters"] else "")

            self.t_entry.delete(0, tk.END)
            self.t_entry.insert(0, model["parameters"].get("T", []) if "T" in model["parameters"] else "")

        elif model.get("type") == "ViscoplasticDesai":
            # Grid setup and entry updates for each parameter
            self.mu_1_entry.delete(0, tk.END)
            self.mu_1_entry.insert(0, model["parameters"].get("mu_1", []) if "mu_1" in model["parameters"] else "")

            self.N_1_entry.delete(0, tk.END)
            self.N_1_entry.insert(0, model["parameters"].get("N_1", []) if "N_1" in model["parameters"] else "")

            self.n_desai_entry.delete(0, tk.END)
            self.n_desai_entry.insert(0, model["parameters"].get("n", []) if "n" in model["parameters"] else "")

            self.a_1_entry.delete(0, tk.END)
            self.a_1_entry.insert(0, model["parameters"].get("a_1", []) if "a_1" in model["parameters"] else "")

            self.eta_desai_entry.delete(0, tk.END)
            self.eta_desai_entry.insert(0, model["parameters"].get("eta", []) if "eta" in model["parameters"] else "")

            self.beta_1_entry.delete(0, tk.END)
            self.beta_1_entry.insert(0, model["parameters"].get("beta_1", []) if "beta_1" in model["parameters"] else "")

            self.beta_entry.delete(0, tk.END)
            self.beta_entry.insert(0, model["parameters"].get("beta", []) if "beta" in model["parameters"] else "")

            self.m_entry.delete(0, tk.END)
            self.m_entry.insert(0, model["parameters"].get("m", []) if "m" in model["parameters"] else "")

            self.gamma_entry.delete(0, tk.END)
            self.gamma_entry.insert(0, model["parameters"].get("gamma", []) if "gamma" in model["parameters"] else "")

            self.alpha_0_entry.delete(0, tk.END)
            self.alpha_0_entry.insert(0, model["parameters"].get("alpha_0", []) if "alpha_0" in model["parameters"] else "")

            self.k_v_entry.delete(0, tk.END)
            self.k_v_entry.insert(0, model["parameters"].get("k_v", []) if "k_v" in model["parameters"] else "")

            self.sigma_t_entry.delete(0, tk.END)
            self.sigma_t_entry.insert(0, model["parameters"].get("sigma_t", []) if "sigma_t" in model["parameters"] else "")

        else:
            self.hide_all_parameters()
        self.loading = False

    def confirm_and_remove(self, category, name):
        response = messagebox.askyesno("Confirm Remove", f"Are you sure you want to remove '{name}' from {category}?")
        if response:
            self.remove_item(category, name)

    def confirm_keep_item_changes(self, category, name):
        response = messagebox.askyesno("Confirm Edit", f"Are you sure you want to Keep Changes '{name}' from {category}?")
        return response

    def remove_item(self, category, name):
        del self.data["constitutive_model"][category][name]
        self.refresh_tree()

    def save_to_file(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[["JSON files", "*.json"]])
        if file_path:
            with open(file_path, "w") as file:
                json.dump(self.data, file, indent=4)

    def load_from_file(self):
        file_path = filedialog.askopenfilename(filetypes=[["JSON files", "*.json"]])
        if file_path:
            with open(file_path, "r") as file:
                self.data = json.load(file)
            self.refresh_tree()

if __name__ == "__main__":
    root = tk.Tk()
    app = JSONConstitutiveApp(root)
    root.mainloop()
