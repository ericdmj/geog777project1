import arcpy
import matplotlib.pyplot as plt
import shapefile  # Library for reading shapefiles
import numpy as np
from arcpy.sa import *
from tkinter import *
from tkinter import ttk, filedialog, font
from tkinter.filedialog import asksaveasfilename
from tkinter import messagebox
from PIL import ImageTk, Image
from mpl_toolkits.basemap import Basemap
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import Normalize, ListedColormap, BoundaryNorm
from matplotlib.figure import Figure


def validate_kentry(value):
    """Validates that the input is a float greater than 1. Shows a messagebox if invalid."""
    try:
        float_value = float(value)
        if float_value > 1:
            return True
        else:
            messagebox.showerror("Invalid Input", "Please enter a number greater than 1.")
            return False
    except ValueError:
        messagebox.showerror("Invalid Input", "Please enter a valid number greater than 1.")
        return False  # Reject non-numeric input


# Create Tkinter window
root = Tk()
root.title("Cancer and Nitrates - A Spatial Analysis")

h1Font = font.Font(family='Garamond', name='appH1Font', size=24, weight='bold')
h2Font = font.Font(family='Garamond', name='appH2Font', size=18, weight='bold')
h3Font = font.Font(family='Garamond', name='appH3Font', size=14, weight='bold')
h4Font = font.Font(family='Helvetica', name='appH4Font', size=10, weight='bold')
bodyFont = font.Font(family='Helvetica', name='appBodyFont', size=10)
radioFont = font.Font(family='Helvetica', name='appRadioFont', size=10)

# Register the validation function
vcmd = root.register(validate_kentry)


# The qualifiedFieldNames environment is used by Copy Features when persisting
# the join field names.
arcpy.env.qualifiedFieldNames = False

# Overwrite existing datasets
arcpy.env.overwriteOutput = True

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

# Set workspace
arcpy.env.workspace = "C:/Users/ericd/Documents/Wisconsin/GEOG777/Project01Python01/shapefiles"

# Initialize K Value variable
kentry_value = StringVar()

# Set local variables
inPointFeatures = "well_nitrate.shp"
zField = "nitr_ran"
cellSize = 0.001
output_raster = "output.tif"

inZoneData = "cancer_tracts.shp"
zoneField = "GEOID10"
inValueRaster = output_raster
outTable = "zonaltable.dbf"

inFeatures = "cancer_tracts"
joinTable = "zonaltable"
joinField = "GEOID10"
outFeature = "zonaljoin"

# Read cancer tract shapefile
cancer_sf = shapefile.Reader("shapefiles/cancer_tracts.shp")
cancer_records = [record.as_dict() for record in cancer_sf.iterRecords()]
canrate_values = np.array([record.get("canrate", 0) for record in cancer_records])

# Define bins for cancer tract colors
cannum_bins = 5
canbins = np.linspace(canrate_values.min(), canrate_values.max(), cannum_bins + 1)
cancer_cmap = plt.get_cmap('Blues')
cancer_norm = BoundaryNorm(canbins, cancer_cmap.N)

# Read well_nitrate shapefile
well_sf = shapefile.Reader("shapefiles/well_nitrate.shp")
well_records = [record.as_dict() for record in well_sf.iterRecords()]
nitr_ran_values = np.array([record.get("nitr_ran", 0) for record in well_records])
well_points = [point.points[0] for point in well_sf.shapes()]  # Extract (x, y)

# Define bins for nitrate values
nitrnum_bins = 5
nitrate_bins = np.linspace(nitr_ran_values.min(), nitr_ran_values.max(), nitrnum_bins + 1)
nitrate_cmap = plt.get_cmap('Oranges')
nitrate_norm = BoundaryNorm(nitrate_bins, nitrate_cmap.N)


# File to save the message
message_file = "glr_output.txt"



# Create Matplotlib figure
fig, ax = plt.subplots(figsize=(10, 6), facecolor="#d9e4d9")

# Create Basemap
map = Basemap(llcrnrlon=-93.13, llcrnrlat=42.23, urcrnrlon=-84.84, urcrnrlat=47.49,
              resolution='i', projection='tmerc', lat_0=44.67, lon_0=-88.85, ax=ax)

# Draw base map features (used in all views)
def draw_basemap():
    map.drawmapboundary(fill_color='#c9ced8')
    map.fillcontinents(color='#799d7a', lake_color='#c9ced8')
    map.drawcoastlines(color='#503b3b')
    map.readshapefile('shapefiles/cancer_tracts', 'cancer_tracts', drawbounds=True)



selected_option = IntVar(value=2)

def threestep():
    try:
        global outbins, glr_message, selected_option
        global canrate_values, mean_nitrate_values, outfeatures_records, outfeatures_values, outfeatures_cmap, outfeatures_norm

        power = kentry_value.get()

        # Path to the output shapefile
        out_features_path = "shapefiles/out_features.shp"

        # Check if the file exists and delete it before creating a new one
        if arcpy.Exists(out_features_path):
            print(f"Deleting existing shapefile: {out_features_path}")
            arcpy.management.Delete(out_features_path)

        # Execute IDW Interpolation
        outIDW = Idw(inPointFeatures, zField, cellSize, power)
        outIDW.save(output_raster)

        messagebox.showinfo("Success", "IDW interpolation completed successfully! Note: the raster file is available in your shapefiles folder as " + output_raster)

        # Execute ZonalStatisticsAsTable
        outZSaT = ZonalStatisticsAsTable(inZoneData, zoneField, inValueRaster,
            outTable, "NODATA", "MEAN", "CURRENT_SLICE", None, None, None, None, None)

        # Join the feature layer to a table
        zonal_joined_table = arcpy.management.AddJoin(inFeatures, joinField, joinTable, joinField)

        # Copy the layer to a new permanent feature class
        result = arcpy.management.CopyFeatures(zonal_joined_table, outFeature)

        messagebox.showinfo("Success yet again", "Your nitrate raster values have been connected to the census tracts! One more step to go.")

        # Perform Generalized Linear Regression
        arcpy.stats.GeneralizedLinearRegression(outFeature, "canrate", "CONTINUOUS", out_features_path, "MEAN")
        glr_message = arcpy.GetMessages()  # Store tool messages

        # Save the message to a file
        with open(message_file, "w", encoding="utf-8") as file:
            file.write(glr_message)

        messagebox.showinfo("Final success!", "And the analysis is done! Note the analysis output buttons on the left are available to you.")

        # After creating the file, read it for visualization
        if arcpy.Exists(out_features_path):
            global outfeatures_records, outfeatures_values, outfeatures_cmap, outfeatures_norm

            outfeatures_sf = shapefile.Reader(out_features_path)
            outfeatures_records = [record.as_dict() for record in outfeatures_sf.iterRecords()]
            outfeatures_values = np.array([record.get("STDRESID", 0) for record in outfeatures_records])

            # Extract Cancer Rate (Y-Axis) and Mean Nitrate Level (X-Axis)
            canrate_values = np.array([record.get("canrate", 0) for record in outfeatures_records])
            mean_nitrate_values = np.array([record.get("MEAN", 0) for record in outfeatures_records])

            outbins = np.array([-3, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3])
            outfeatures_cmap = plt.get_cmap('PRGn')
            outfeatures_norm = BoundaryNorm(outbins, outfeatures_cmap.N)
            #print(f"Threestep Bins: {outbins}")

        # Automatically refresh the map if Standard Residuals is selected
        if selected_option.get() == 3:
            update_map()

        # Enable the Standard Residuals radio button after successful execution
        std_res_radio.config(state="normal")
        open_button.config(state="normal")

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred:\n{e}")






# Track colorbars globally to remove them before adding new ones
colorbars = []

# Track title and description labels globally
map_title = None
map_description = None

def update_map():
    """Update the displayed map based on the selected radio button."""
    global colorbars, ax, map_title, map_description, outbins, glr_message, save_map_file

    # Clear the figure to prevent shrinking
    fig.clear()
    ax = fig.add_subplot(111)

    # Recreate Basemap to maintain correct scaling
    map = Basemap(llcrnrlon=-93.13, llcrnrlat=42.23, urcrnrlon=-84.84, urcrnrlat=47.49,
                  resolution='i', projection='tmerc', lat_0=44.67, lon_0=-88.85, ax=ax)

    # Remove previous colorbars safely
    for cb in colorbars:
        if cb.ax in fig.axes:  # Only remove if still in the figure
            cb.remove()
    colorbars.clear()

    # Remove previous title and description if they exist
    if map_title is not None:
        map_title.destroy()
    if map_description is not None:
        map_description.destroy()

    # Draw base map features
    map.drawmapboundary(fill_color='#c9ced8')
    map.fillcontinents(color='#799d7a', lake_color='#c9ced8')
    map.drawcoastlines(color='#503b3b')
    map.readshapefile('shapefiles/cancer_tracts', 'cancer_tracts', drawbounds=True)

    selected_option = map_option.get()

    if selected_option == 0:  # Cancer Tracts Choropleth Only
        for record, shape, canrate in zip(cancer_records, map.cancer_tracts, canrate_values):
            bin_index = np.digitize(canrate, canbins, right=True) - 1
            color = cancer_cmap(bin_index / (cannum_bins - 1))
            x, y = zip(*shape)
            ax.fill(x, y, color=color, edgecolor="#503b3b")

        map_title = ttk.Label(right_frame, text='Cancer Rates', style='R.TLabel', font=h3Font, foreground="#47466a")
        map_description = ttk.Label(right_frame, text='Cancer rate in Wisconsin by census tract', style='R.TLabel', font=h4Font, foreground="#4f3534")
        save_map_file="cancer_tracts.png"

        sm_cancer = plt.cm.ScalarMappable(cmap=cancer_cmap, norm=cancer_norm)
        sm_cancer.set_array([])
        cbar_cancer = plt.colorbar(sm_cancer, orientation="vertical", label="Cancer Rate", ax=ax, shrink=0.8)
        colorbars.append(cbar_cancer)

    elif selected_option == 1:  # Well Nitrate Points Only
        for well, nitr_ran in zip(well_points, nitr_ran_values):
            x, y = map(well[0], well[1])
            bin_index = np.digitize(nitr_ran, nitrate_bins, right=True) - 1
            color = nitrate_cmap(bin_index / (nitrnum_bins - 1))
            ax.scatter(x, y, color=color, edgecolors='#503b3b', s=15, alpha=0.8)

        map_title = ttk.Label(right_frame, text='Nitrate Levels', style='R.TLabel', font=h3Font, foreground="#47466a")
        map_description = ttk.Label(right_frame, text='Collected from test wells throughout the state', style='R.TLabel', font=h4Font, foreground="#4f3534")
        save_map_file="well_points.png"

        sm_nitrate = plt.cm.ScalarMappable(cmap=nitrate_cmap, norm=nitrate_norm)
        sm_nitrate.set_array([])
        cbar_nitrate = plt.colorbar(sm_nitrate, orientation="vertical", label="Nitrate Level", ax=ax, shrink=0.8)
        colorbars.append(cbar_nitrate)


    elif selected_option == 2:  # Full Map (Choropleth + Nitrate Points)
        for record, shape, canrate in zip(cancer_records, map.cancer_tracts, canrate_values):
            bin_index = np.digitize(canrate, canbins, right=True) - 1
            color = cancer_cmap(bin_index / (cannum_bins - 1))
            x, y = zip(*shape)
            ax.fill(x, y, color=color, edgecolor="#503b3b")

        for well, nitr_ran in zip(well_points, nitr_ran_values):
            x, y = map(well[0], well[1])
            bin_index = np.digitize(nitr_ran, nitrate_bins, right=True) - 1
            color = nitrate_cmap(bin_index / (nitrnum_bins - 1))
            ax.scatter(x, y, color=color, edgecolors='#503b3b', s=15, alpha=0.8)

        map_title = ttk.Label(right_frame, text='Combined Nitrates and Cancer Rates', style='R.TLabel', font=h3Font, foreground="#47466a")
        map_description = ttk.Label(right_frame, text='Nitrate wells by test location superimposed on cancer rates by census tract', style='R.TLabel', font=h4Font, foreground="#4f3534")
        save_map_file="wells_cancer.png"

        sm_cancer = plt.cm.ScalarMappable(cmap=cancer_cmap, norm=cancer_norm)
        sm_cancer.set_array([])
        cbar_cancer = plt.colorbar(sm_cancer, orientation="vertical", label="Cancer Rate", ax=ax, shrink=0.8)
        colorbars.append(cbar_cancer)

        sm_nitrate = plt.cm.ScalarMappable(cmap=nitrate_cmap, norm=nitrate_norm)
        sm_nitrate.set_array([])
        cbar_nitrate = plt.colorbar(sm_nitrate, orientation="vertical", label="Nitrate Level", ax=ax, shrink=0.8)
        colorbars.append(cbar_nitrate)


    elif selected_option == 3:  # Standard Residuals Choropleth Only
        if not arcpy.Exists("shapefiles/out_features.shp"):
            messagebox.showwarning("Missing Data", "Standard Residual data is not available yet. Run analysis first.")
            return

        for record, shape, STDRESID in zip(outfeatures_records, map.cancer_tracts, outfeatures_values):
            bin_index = np.digitize(STDRESID, outbins, right=True) - 1
            #print(f"SOURCE_ID_1: {record.get('SOURCE_ID')} | STDRESID: {STDRESID} | Bin Index Before Clipping: {bin_index}")
            color = outfeatures_cmap(bin_index / (len(outbins) - 2))
            x, y = zip(*shape)
            ax.fill(x, y, color=color, edgecolor="#503b3b")
            #print(f"SOURCE_ID_2: {record.get('SOURCE_ID')} | STDRESID: {STDRESID} | Bin Index After Clipping: {bin_index} | Expected Color: {color}")

        map_title = ttk.Label(right_frame, text='Standard Residuals', style='R.TLabel', font=h3Font, foreground="#47466a")
        map_description = ttk.Label(right_frame, text='After the regression analysis you just performed, the residuals are the differences between observed values and estimated values. Remember, extreme residuals indicate poor model fit.', wraplength=650, style='R.TLabel', font=h4Font, foreground="#4f3534")
        save_map_file="output.png"

        sm_outfeatures = plt.cm.ScalarMappable(cmap=outfeatures_cmap, norm=outfeatures_norm)
        sm_outfeatures.set_array([])
        cbar_outfeatures = plt.colorbar(sm_outfeatures, orientation="vertical", label="Standardized Residuals", ax=ax, shrink=0.8, ticks=outbins[1:-1])
        cbar_outfeatures.ax.set_yticklabels(["-2.5", "-1.5", "-0.5", "0.5", "1.5", "2.5"])
        colorbars.append(cbar_outfeatures)

    map_title.grid(column=0, row=0, columnspan=3, ipady=5, ipadx=5, pady=(20,0))
    map_description.grid(column=0, row=1, columnspan=3, ipadx=5)

    # Refresh the canvas
    canvas.draw()


# Function to open a new window and display the message
def open_message_window():
    global canrate_values, mean_nitrate_values

    if 'mean_nitrate_values' not in globals() or 'canrate_values' not in globals():
        messagebox.showerror("Error", "Data not found. Run the analysis first!")
        return  # Exit if data is missing

    new_window = Toplevel(root)
    new_window.title("GLR Tool Output")

    # Create Notebook (Tabbed Interface)
    notebook = ttk.Notebook(new_window)
    notebook.grid(row=0, column=0, sticky="nsew")

    # Ensure the window resizes properly
    new_window.columnconfigure(0, weight=1)
    new_window.rowconfigure(0, weight=1)

    # === First Tab: GLR Output Message ===
    frame_message = ttk.Frame(notebook)
    notebook.add(frame_message, text="GLR Output Message")

    frame_message.columnconfigure(0, weight=1)
    frame_message.rowconfigure(0, weight=1)

    text_box = Text(frame_message, height=42, width=159, wrap="word")
    text_box.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

    try:
        with open(message_file, "r", encoding="utf-8") as file:
            saved_message = file.read()
        text_box.insert("1.0", saved_message)
    except FileNotFoundError:
        text_box.insert("1.0", "No saved message found.")

    text_box.config(state="disabled")  # Make read-only

    # === "Save As" Button for GLR Output Message ===
    def save_message():
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")],
            initialfile="glr_output.txt",
            title="Save GLR Output As"
        )
        if file_path:
            with open(file_path, "w", encoding="utf-8") as file:
                file.write(saved_message)
            messagebox.showinfo("Success", f"GLR Output saved at:\n{file_path}")

    save_message_button = ttk.Button(frame_message, text="Save As", command=save_message)
    save_message_button.grid(row=1, column=0, pady=5)

    # === Second Tab: GLR Graph (Scatterplot) ===
    frame_graph = ttk.Frame(notebook)
    notebook.add(frame_graph, text="GLR Graph")

    frame_graph.columnconfigure(0, weight=1)
    frame_graph.rowconfigure(0, weight=1)

    fig = Figure(figsize=(8, 6))
    ax1 = fig.add_subplot(111)

    ax1.scatter(mean_nitrate_values, canrate_values, alpha=0.6, edgecolor="black", label="Data Points")
    ax1.set_title("Cancer Rate vs. Nitrate Level")
    ax1.set_xlabel("Mean Nitrate Level (MEAN)")
    ax1.set_ylabel("Cancer Rate (canrate)")

    # Compute R² if there are enough data points
    if len(mean_nitrate_values) > 1:
        # Fit a linear regression model
        m, b = np.polyfit(mean_nitrate_values, canrate_values, 1)

        # Compute predicted values
        y_pred = m * mean_nitrate_values + b

        # Compute R²
        ss_total = np.sum((canrate_values - np.mean(canrate_values))**2)
        ss_residual = np.sum((canrate_values - y_pred)**2)
        r_squared = 1 - (ss_residual / ss_total)

        # Plot the regression line
        ax1.plot(mean_nitrate_values, y_pred, color="red", linestyle="--", label=f"Trendline (R²={r_squared:.3f})")

        # Add R² as a text label inside the plot
        ax1.text(
            0.05, 0.90, f"R² = {r_squared:.3f}",
            transform=ax1.transAxes, fontsize=12,
            verticalalignment="top", bbox=dict(facecolor="white", alpha=0.6)
        )

    ax1.legend()

    canvas = FigureCanvasTkAgg(fig, master=frame_graph)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

    # === "Save As" Button for Scatterplot ===
    def save_scatterplot():
        file_path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG Files", "*.png"), ("All Files", "*.*")],
            initialfile="scatterplot.png",
            title="Save Scatterplot As"
        )
        if file_path:
            fig.savefig(file_path, dpi=300, bbox_inches="tight")
            messagebox.showinfo("Success", f"Scatterplot saved at:\n{file_path}")

    save_scatter_button = ttk.Button(frame_graph, text="Save As", command=save_scatterplot)
    save_scatter_button.grid(row=1, column=0, pady=5)




# Style the app
s = ttk.Style()
s.configure('Main.TFrame', background='#bccebc', borderwidth=15, relief='ridge')
s.configure("TButton", padding=3, relief="flat", foreground="#AA0000", background="#bccebc", font=bodyFont)
s.configure('TSeparator', background="#bccebc")
s.configure('TLabel', background="#bccebc")
s.configure('R.TLabel', background="#d9e4d9")
s.configure('L.TRadiobutton', font=radioFont, foreground="#4f3534", background="#d9e4d9", padding=(10,10,10,10))
s.configure('K.TFrame', background='#bccebc', borderwidth=2, relief='ridge')

# Tkinter Grid Layout
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

# Create frames for the map
main_frame = ttk.Frame(root, padding=(10, 10, 10, 10), style="Main.TFrame")
main_frame.grid(row=0, column=0, sticky="nsew")

top_frame = Frame(main_frame, background="#bccebc")
top_frame.columnconfigure(0, weight=1)
left_frame = Frame(main_frame, background="#bccebc")
right_frame = Frame(main_frame, background="#d9e4d9")

top_frame.grid(column=0, row=0, columnspan=2, sticky=(N,S,E,W), padx=20)
left_frame.grid(column=0, row=1, sticky=(N,S,E,W), padx=20)
right_frame.grid(column=1, row=1, sticky=(N,S,E,W), padx=20)

topframe_pair = Frame(top_frame, background="#bccebc")
topframe_pair.grid(column=0, row=0, sticky=(W))

# Create the different widgets; widgets in the lefthand frame
toptitlelogo = ImageTk.PhotoImage(Image.open('nitrate2.png'))
logo_label = ttk.Label(topframe_pair, image=toptitlelogo)
toptitletext = ttk.Label(topframe_pair, text='Cancer and Nitrates - A Spatial Analysis', font=h1Font, foreground="#4f3534")
hbar = ttk.Separator(top_frame, orient='horizontal')
maintitletext = ttk.Label(left_frame, text='Explore', font=h2Font, foreground="#47466a")
mainexplanatory = ttk.Label(left_frame, text="A. Check out the Base Shapefiles", font=h4Font, foreground="#4f3534", justify="center")

kheader = ttk.Label(left_frame, text='B. Inverse Distance Weighting Test', font=h4Font, foreground="#4f3534", justify="center")
kentry = ttk.Label(left_frame, text="Let's explore the relationship between nitrate and cancer through spatial analysis. What distance exponent 'k' would you like to use?", wraplength=175, font=bodyFont, foreground="#4f3534")
kentry_box = ttk.Entry(left_frame, textvariable=kentry_value, width=6, validate="focusout", validatecommand=(vcmd, "%P"), justify="center", font=h3Font, foreground="#47466a")
ksubmit = ttk.Button(left_frame, text='Submit', default='active', command=threestep)


# Radio buttons frame
map_option = IntVar(value=2)  # Default to full map
radio_frame = Frame(left_frame, background="#bccebc")
radio_frame.grid(row=2, column=0, pady=10, padx=5)
rb1 = ttk.Radiobutton(radio_frame, text="Full Map (Both)", variable=map_option, value=2, command=update_map, style="L.TRadiobutton")
rb1.grid(row=0, column=0, sticky="ew", pady=5)
rb2 = ttk.Radiobutton(radio_frame, text="Cancer Tracts Only", variable=map_option, value=0, command=update_map, style="L.TRadiobutton")
rb2.grid(row=1, column=0, sticky="ew", pady=5)
rb3 = ttk.Radiobutton(radio_frame, text="Well Nitrate Points Only", variable=map_option, value=1, command=update_map, style="L.TRadiobutton")
rb3.grid(row=2, column=0, sticky="ew", pady=5)


c_header = ttk.Label(left_frame, text='C. Explore the Results through\nLinear Regression', font=h4Font, foreground="#4f3534", justify="center")
c_entry = ttk.Label(left_frame, text="After the analysis is run, you can view the Standard Residuals map and view the output report + graph", wraplength=175, font=bodyFont, foreground="#4f3534")


std_res_radio = ttk.Radiobutton(left_frame, text="Standard Residuals", variable=map_option, value=3, command=update_map, style="L.TRadiobutton")
std_res_radio.grid(row=9, column=0, sticky="ew", pady=(10,5), padx=22)

# Start with the Standard Residuals button disabled
std_res_radio.config(state="disabled")

# Create a button to open the new window
open_button = ttk.Button(left_frame, text="View GLR Output", command=open_message_window)
open_button.grid(row=10, column=0, padx=20, pady=10)
open_button.config(state="disabled")

# Grid the widgets
logo_label.grid(column=0, row=0,padx=3, pady=3)
toptitletext.grid(column=1, row=0, sticky=W, padx=10, pady=10)
hbar.grid(column=0, row=1, sticky='ew', columnspan=2, pady=10)
maintitletext.grid(column=0, row=0, pady=10)
mainexplanatory.grid(column=0, row=1, pady=10)

kheader.grid(column=0, row=3, pady=10)
kentry.grid(column=0, row=4, pady=(10, 0))
kentry_box.grid(column=0, row=5, pady=10)
ksubmit.grid(column=0, row=6, pady=(5,10))

c_header.grid(column=0, row=7, pady=10)
c_entry.grid(column=0, row=8, pady=(10, 0))



# Embed Matplotlib figure in Tkinter
canvas = FigureCanvasTkAgg(fig, master=right_frame)
canvas.draw()
canvas.get_tk_widget().grid(row=3, column=0, sticky="nsew")


# Action buttons
action_frame = Frame(right_frame, background="#d9e4d9")
action_frame.grid(row=6, column=0, pady=10)



# Function to save the message to a different file
def save_map():
    """Save the currently displayed map as a PNG file using a Save As dialog with a suggested filename."""
    # Suggested filename based on the currently selected map option
    map_names = {
        0: "cancer_tracts.png",
        1: "well_nitrate_points.png",
        2: "full_map.png",
        3: "standard_residuals.png"
    }
    suggested_filename = map_names.get(map_option.get(), "map.png")  # Default to "map.png" if unknown

    file_path = filedialog.asksaveasfilename(
        defaultextension=".png",
        filetypes=[("PNG Files", "*.png"), ("All Files", "*.*")],
        initialfile=suggested_filename,  # Suggested file name
        title="Save Map As"
    )

    if file_path:  # Proceed only if a file path is selected
        canvas.figure.savefig(file_path, dpi=300, bbox_inches='tight')
        messagebox.showinfo("Success", f"Map saved successfully at:\n{file_path}")



# Add a "Save As" button
save_button = ttk.Button(action_frame, text="Save This Map", command=save_map)
save_button.grid(column=0, row=1, pady=5)



# Define closing function to release resources
def on_closing():
    plt.close(fig)  # Close the Matplotlib figure
    root.quit()     # Quit Tkinter event loop
    root.destroy()  # Destroy the Tkinter window

# Bind the closing function to the window close button
root.protocol("WM_DELETE_WINDOW", on_closing)


update_map()  # Initialize the map
root.mainloop()
