import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox
import os
import io
import contextlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class PythonScriptRunner(tk.Frame):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.pack(fill="both", expand=True)
        self.script_path = None

        # 2 frame
        self.top_frame = tk.Frame(self, height=560)  # 70% of 800
        self.top_frame.pack(side="top", fill="both", expand=False)
        self.top_frame.pack_propagate(False) 

        self.bottom_frame = tk.Frame(self, height=290)  # 30% of 800
        self.bottom_frame.pack(side="top", fill="both", expand=False)
        self.bottom_frame.pack_propagate(False)

        # top frame

        # file selection
        choose_frame = tk.Frame(self.top_frame)
        choose_frame.pack(fill="x", padx=10, pady=5)

        self.choose_button = tk.Button(choose_frame, text="Choose Python Script", command=self.choose_script)
        self.choose_button.pack(side="left")

        self.script_label = tk.Label(choose_frame, text="No file selected", anchor="w")
        self.script_label.pack(side="left", fill="x", expand=True, padx=10)

        self.run_button = tk.Button(choose_frame, text="Run", command=self.run_script)
        self.run_button.pack(side="right")

        # code
        code_label = tk.Label(self.top_frame, text="Script Content:")
        code_label.pack(anchor="w", padx=10)

        self.code_text = scrolledtext.ScrolledText(self.top_frame)
        self.code_text.pack(fill="both", expand=True, padx=10, pady=5)

        # bottom_frame

        output_label = tk.Label(self.bottom_frame, text="Output:")
        output_label.pack(anchor="w", padx=10)

        self.output_text = scrolledtext.ScrolledText(self.bottom_frame)
        self.output_text.pack(fill="both", expand=True, padx=10, pady=5)

        # diagram
        self.plot_frame = tk.Frame(self)
        self.plot_frame.pack(padx=10, pady=5, fill="both", expand=True)

    def choose_script(self):
        file_path = filedialog.askopenfilename(filetypes=[("Python files", "*.py")])
        if file_path:
            self.script_path = file_path
            self.script_label.config(text=os.path.basename(file_path))
            with open(file_path, 'r', encoding='utf-8') as file:
                code = file.read()
                self.code_text.delete(1.0, tk.END)
                self.code_text.insert(tk.END, code)

    def run_script(self):
        if not self.script_path:
            messagebox.showwarning("No File", "Please choose a Python script first.")
            return

        self.output_text.delete(1.0, tk.END)
        self.clear_plots()

        try:
            with open(self.script_path, 'r', encoding='utf-8') as file:
                script_code = file.read()

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                plt.close('all')
                exec(script_code, {'__name__': '__main__'})

            output = buffer.getvalue()
            self.output_text.insert(tk.END, output)

            self.draw_matplotlib_figures()

        except Exception as e:
            self.output_text.insert(tk.END, f"Error: {str(e)}")

    def clear_plots(self):
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

    def draw_matplotlib_figures(self):
        figures = [plt.figure(n) for n in plt.get_fignums()]
        for fig in figures:
            canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(side="top", fill="both", expand=True)


if __name__ == "__main__":
    root = tk.Tk()
    root.title("Standalone Python Script Runner")
    root.geometry("1000x800")
    runner = PythonScriptRunner(master=root)
    root.mainloop()
