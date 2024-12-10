import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import math
import numpy as np
import csv
import ast
import subprocess
import shutil
import argparse

import os
os.environ['PATH'] += os.pathsep + r'C:\Program Files\Graphviz\bin'
import graphviz

from PIL import Image, ImageDraw, ImageFont

parser = argparse.ArgumentParser(description="Process some integers.")
parser.add_argument('rep', type=int, help='An integer specifying the repetition number')
args = parser.parse_args()
rep = args.rep


arrowSize = .4
nodeSize = 1.25


renderFragMaps = True
renderFullDataFlow = True
renderStepsDataFlow = True
renderState2State = True
renderUpdates = 240#27*12 # 240 # 2400 per lifetime!
print("renderUpdates = ",renderUpdates)

features = ['direction','front_open']
bg_color = 'white'
###########################################################################
## LOAD CONFIG
###########################################################################

def load_csv_with_lists(file): # file includes path
    data = []
    with open(file, mode='r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Convert string representations of lists to actual lists
            for key, value in row.items():
                try:
                    # Use ast.literal_eval to safely evaluate strings containing Python literals
                    row[key] = ast.literal_eval(value) if value.startswith('[') else value
                except (ValueError, SyntaxError):
                    # If parsing fails, keep the original string
                    pass
            data.append(row)
    return data

def convert_data(data):
    result = {}
    for item in data:
        key = item['name']
        value = item['values']
        # Try to split by comma and convert to integer list if applicable
        if ',' in value:
            try:
                value = [int(x) for x in value.split(',')]
            except ValueError:
                pass  # If conversion fails, keep the original value
        elif value.isdigit():  # If the value is a number, convert to int
            value = int(value)
        
        result[key] = value
    return result

###########################################################################
## GRAPHVIZ TOOLS
###########################################################################


def rgb_to_hex(rgb):
    """
    Convert an RGB list [R, G, B] where each value is between 0 and 1
    into a hexadecimal string '#XXXXXX'.
    """
    return '#{:02X}{:02X}{:02X}'.format(
        int(min(1.0,rgb[0]) * 255), 
        int(min(1.0,rgb[1]) * 255), 
        int(min(1.0,rgb[2]) * 255)
    )


def convert_dot_txt_to_png(inFile,outFile):
    # Read the DOT content from the text file
    with open(inFile, 'r') as file:
        dot_content = file.read()
    
    # Create a graph from the DOT content
    graph = graphviz.Source(dot_content)
    
    # Render and save the graph as a PNG file
    graph.render(outFile[:-4], format='png', cleanup=True)


def create_graphviz_rich(neighbor_dict,values_in_order,ids_in_order,config,features,connectionWeighting = 5):
    for f_id in range(len(features)):
        feature = features[f_id]
        print('    feature:',feature)
        graph_lines = ["digraph G {"]  # Start of the Graphviz DOT format
        graph_lines.append('bgcolor="black";  // Sets the background color to black')
        graph_lines.append('node [shape=ellipse, fixedsize=true, width='+str(nodeSize) + ', height='+str(.75 * nodeSize)+'];')
        graph_lines.append('node [label=""];')
        count = 0
        for key, neighbors in neighbor_dict.items():
    
            if isinstance(neighbors, int):
                neighbors = [neighbors]  # Convert x into a list containing x
            if neighbors[-1] != key:
                neighbors += [key]  # Convert x into a list containing x
    
            if isinstance(config['T0_cells'], int):
                config['T0_cells'] = [config['T0_cells']]  # Convert x into a list containing x
            if isinstance(config['T1_cells'], int):
                config['T1_cells'] = [config['T1_cells']]  # Convert x into a list containing x
    
            
            cc = [values_in_order[f_id][ids_in_order.index(key)] * 10,
                  values_in_order[f_id][ids_in_order.index(key)] * 10,
                  values_in_order[f_id][ids_in_order.index(key)] * 10]
            fontColor = "white"
            if sum(cc) > .33:
                fontColor = "black"
            if key in config['T0_cells']:
                graph_lines.append(str(key) + ' [style = "filled",fillcolor = "' + rgb_to_hex(cc) + '",fontcolor  = "'+fontColor+'", color="blue", penwidth=20]')
            elif key in config['T1_cells']:
                graph_lines.append(str(key) + ' [style = "filled",fillcolor = "' + rgb_to_hex(cc) + '",fontcolor  = "'+fontColor+'", color="orange", penwidth=20]')
            else:   
                graph_lines.append(str(key) + ' [style = "filled",fillcolor = "' + rgb_to_hex(cc) + '",fontcolor  = "'+fontColor+'", color="white", penwidth=3]')
            for i, neighbor in enumerate(neighbors):
                # Add edge from key to neighbor
                # load this cells predictions!

                this_fragMap=load_frag_map(inPath+'vis_output/cellPredictions/cellPrediction_'+str(key)+'.py')

                #if step == None:
                #    this_fragMap=load_frag_map(path+'vis_output/cellPredictions/cellPrediction_'+str(key)+'.py')
                #else:
                #    this_fragMap=load_frag_map(path+'vis_output/cellPredictions/cellPrediction_'+str(key)+'_step_'+str(step)+'.py')
                #    #print("loading: ", key, step)
    
                
                #graph_lines.append(f"    {key} -> {neighbor};")
                penWidth = this_fragMap[0][i] * connectionWeighting
                if penWidth < .2:
                    penWidth = 0
                else:
                    penWidth = 3 + 2 * penWidth
                if penWidth > 100:
                    penWidth = 100

                #graph_lines.append(f'    {neighbor} -> {key} [penwidth={penWidth}, color="white", arrowsize = '+str(penWidth*arrowSize)+'];')
                graph_lines.append(f'    {neighbor} -> {key} [penwidth={penWidth}, color="white", arrowsize = '+str(.5+.125*penWidth*arrowSize)+'];')
                #print(graph_lines[-1])
                # Add edge from neighbor to key
                #graph_lines.append(f"    {neighbor} -> {key};")
            ## # add self
            ## penWidth = min(this_fragMap[0][len(neighbors)] * 5,100)
            ## if penWidth < .001:
            ##     penWidth = 0
            ## graph_lines.append(f'    {key} -> {key} [penwidth={penWidth}, color="white"];')
    
        count += 1
        graph_lines.append("}")  # End of the Graphviz DOT format
        graphviz_text = "\n".join(graph_lines)

        #    file.write(graphviz_text)
        tempFile = tempPath+'temp.txt'               # name of temp file with graphviz code
        outFile = outPath+'FLOW_FULL_'+feature+'.png' # name of file with rendered graph
        with open(tempFile, 'w') as file:
            file.write(graphviz_text)
        convert_dot_txt_to_png(tempFile,outFile)
        
        # Load the image
        img = Image.open(outFile)

        # Specify the height of the additional space for text
        additional_height = 400  # Adjust as needed for the amount of text
        
        # Create a new image with extra space at the bottom
        new_img = Image.new('RGB', (img.width, img.height + additional_height), (0, 0, 0))
        new_img.paste(img, (0, 0))

        # Prepare to draw on the image
        draw = ImageDraw.Draw(new_img)
        
        # Optional: Load a font (or use a default one)
        font = ImageFont.truetype("arial.ttf", size=50)  # You can specify your own font
        #font = ImageFont.load_default()
        
        # Define the text and its position
        text = "Feature: " + feature +"     steps shown as green -> blue -> red (early to late)\n"
        text += "sensor input: " + str(config['T0_cells'])+"\n\n"
        position = (10, new_img.height - 300)  # Adjust position as needed
        draw.text(position, text, font=font, size=50, fill="white")

        text = "\n"
        text += "sensor wall bump: " + str(config['T0_cells'][0])+"\n"
        text += "sensor wall front: " + str(config['T0_cells'][1])+"\n"
        text += "sensor direction A: " + str(config['T0_cells'][2])+"\n"
        text += "sensor direction B: " + str(config['T0_cells'][3])+"\n"
        position = (new_img.width * .25, new_img.height - 300)  # Adjust position as needed
        draw.text(position, text, font=font, size=50, fill="white")
        
        # Save or display the modified image
        new_img.save(outFile)
        new_img.close()


def create_graphviz_rich_steps(neighbor_dict,values_in_order,ids_in_order,config,features,timeColorConnections=False,connectionWeighting = 5):
    
    neighbor_ids,ids_in_order = load_neighbor_ids(inPath+'vis_output/neighbor_ids.txt')

    for f_id in range(len(features)):
        feature = features[f_id]
        print('    feature:',feature,end=' ')
        graph_lines = ["digraph G {"]  # Start of the Graphviz DOT format
        graph_lines.append('bgcolor="black";  // Sets the background color to black')
        graph_lines.append('node [shape=ellipse, fixedsize=true, width='+str(nodeSize) + ', height='+str(.75 * nodeSize)+'];')
        graph_lines.append('node [label=""];')

        count = 0
        for key, neighbors in neighbor_dict.items():
            print('.',end='')
    
            if isinstance(neighbors, int):
                neighbors = [neighbors]  # Convert x into a list containing x
            if neighbors[-1] != key:
                neighbors += [key]  # Convert x into a list containing x
    
            if isinstance(config['T0_cells'], int):
                config['T0_cells'] = [config['T0_cells']]  # Convert x into a list containing x
            if isinstance(config['T1_cells'], int):
                config['T1_cells'] = [config['T1_cells']]  # Convert x into a list containing x

            cc = [0,0,0]
            for i in range(12):
                cc[0] += values_in_order[i][f_id][ids_in_order.index(key)] * max(0.0,float(i-7))
                cc[1] += values_in_order[i][f_id][ids_in_order.index(key)] * max(0.0,float(4-i))
                cc[2] += values_in_order[i][f_id][ids_in_order.index(key)] * max(0.0,((12.0 - 2.182 * abs(float(i) - 5.5)) - 3) * 4/3)
            cc[0]*=.2;
            cc[1]*=.2;
            cc[2]*=.1;
            
            fontColor = "white"
            if sum(cc) > .50:
                fontColor = "black"
            if key in config['T0_cells']:
                graph_lines.append(str(key) + ' [style = "filled",fillcolor = "' + rgb_to_hex(cc) + '",fontcolor  = "'+fontColor+'", color="blue", penwidth=20]')
            elif key in config['T1_cells']:
                graph_lines.append(str(key) + ' [style = "filled",fillcolor = "' + rgb_to_hex(cc) + '",fontcolor  = "'+fontColor+'", color="orange", penwidth=20]')
            else:   
                graph_lines.append(str(key) + ' [style = "filled",fillcolor = "' + rgb_to_hex(cc) + '",fontcolor  = "'+fontColor+'", color="#A0A0A0", penwidth=3]')
            for i, neighbor in enumerate(neighbors):
                # Add edge from key to neighbor
                # load this cells predictions!
                
                this_fragMap=load_frag_map(inPath+'vis_output/cellPredictions/cellPrediction_'+str(key)+'.py')
                
                #if step == None:
                #    this_fragMap=load_frag_map(path+'vis_output/cellPredictions/cellPrediction_'+str(key)+'.py')
                #else:
                #    this_fragMap=load_frag_map(path+'vis_output/cellPredictions/cellPrediction_'+str(key)+'_step_'+str(step)+'.py')
                #    #print("loading: ", key, step)

                linkCC = [0,0,0]
                for s in range(12):

                    link_fragMap = load_frag_map(inPath+'vis_output/cellPredictions/cellPrediction_'+str(key)+'_step_'+str(s)+'.py')
                    #print(link_fragMap)
                    #print("s,key",s,key)
                    #print("len(ids_in_order)",len(ids_in_order))
                    #print("len(link_fragMap[0])",len(link_fragMap[0]))
                    #print(ids_in_order.index(key),ids_in_order.index(key))
                    
                    linkCC[0] += (link_fragMap[0][i] * max(0.0,float(s-7)))
                    linkCC[1] += (link_fragMap[0][i] * max(0.0,float(4-s)))
                    linkCC[2] += (link_fragMap[0][i] * max(0.0,((12.0 - 2.182 * abs(float(s) - 5.5)) - 3) * 4/3))
                    
                linkCC[0]*=.2;
                linkCC[1]*=.2;
                linkCC[2]*=.1;

                linkCC[0]+=.1;
                linkCC[1]+=.1;
                linkCC[2]+=.1;
                
                #graph_lines.append(f"    {key} -> {neighbor};")
                penWidth = this_fragMap[0][i] * connectionWeighting
                if penWidth < .2:
                    penWidth = 0
                else:
                    penWidth = 3 + 2 * penWidth
                if penWidth > 100:
                    penWidth = 100
    
                if timeColorConnections:
                    if sum(linkCC) > .31:
                        graph_lines.append(f'    {neighbor} -> {key} [penwidth={penWidth}, color="' + rgb_to_hex(linkCC) + '", arrowsize = '+str(.5+.125*penWidth*arrowSize)+'];')
                    else:
                        graph_lines.append(f'    {neighbor} -> {key} [penwidth=0, color="' + rgb_to_hex(linkCC) + '", arrowsize = '+str(.5+.125*penWidth*arrowSize)+'];')
                else:
                    graph_lines.append(f'    {neighbor} -> {key} [penwidth={penWidth}, color="white", arrowsize = '+str(.5+.125*penWidth*arrowSize)+'];')
                #print(graph_lines[-1])
                # Add edge from neighbor to key
                #graph_lines.append(f"    {neighbor} -> {key};")
            ## # add self
            ## penWidth = min(this_fragMap[0][len(neighbors)] * 5,100)
            ## if penWidth < .001:
            ##     penWidth = 0
            ## graph_lines.append(f'    {key} -> {key} [penwidth={penWidth}, color="white"];')

        print('  generating plot')
        count += 1
        graph_lines.append("}")  # End of the Graphviz DOT format
        graphviz_text = "\n".join(graph_lines)

        #    file.write(graphviz_text)
        tempFile = tempPath+'temp.txt'               # name of temp file with graphviz code
        if timeColorConnections:
            outFile = outPath+'FLOW_STEP_COLOR_'+feature+'.png' # name of file with rendered graph
        else:
            outFile = outPath+'FLOW_STEP_'+feature+'.png' # name of file with rendered graph

        with open(tempFile, 'w') as file:
            file.write(graphviz_text)
        convert_dot_txt_to_png(tempFile,outFile)
        # Load the image
        
        img = Image.open(outFile)

        # Specify the height of the additional space for text
        additional_height = 400  # Adjust as needed for the amount of text
        
        # Create a new image with extra space at the bottom
        new_img = Image.new('RGB', (img.width, img.height + additional_height), (0, 0, 0))
        new_img.paste(img, (0, 0))

        # Prepare to draw on the image
        draw = ImageDraw.Draw(new_img)
        
        # Optional: Load a font (or use a default one)
        font = ImageFont.truetype("arial.ttf", size=50)  # You can specify your own font
        #font = ImageFont.load_default()
        
        # Define the text and its position
        text = "Feature: " + feature +"     steps shown as green -> blue -> red (early to late)\n"
        text += "sensor input: " + str(config['T0_cells'])+"\n\n"
        position = (10, new_img.height - 300)  # Adjust position as needed
        draw.text(position, text, font=font, size=50, fill="white")

        text = "\n"
        text += "sensor wall bump: " + str(config['T0_cells'][0])+"\n"
        text += "sensor wall front: " + str(config['T0_cells'][1])+"\n"
        text += "sensor direction A: " + str(config['T0_cells'][2])+"\n"
        text += "sensor direction B: " + str(config['T0_cells'][3])+"\n"
        position = (new_img.width * .25, new_img.height - 300)  # Adjust position as needed
        draw.text(position, text, font=font, size=50, fill="white")
        
        # Save or display the modified image
        new_img.save(outFile)
        new_img.close()


def create_graphviz_brainState(neighbor_dict,values_in_order,states_in_order,ids_in_order,config,connectionWeighting = 5):
    graph_lines = ["digraph G {"]  # Start of the Graphviz DOT format
    graph_lines.append('node [shape=ellipse, fixedsize=true, width='+str(nodeSize) + ', height='+str(.75 * nodeSize)+'];')
    #graph_lines.append('node [label=""];')
    graph_lines.append('bgcolor="'+bg_color+'";  // Sets the background color')

    count = 0
    for key, neighbors in neighbor_dict.items():

        if isinstance(neighbors, int):
            neighbors = [neighbors]  # Convert x into a list containing x
        if neighbors[-1] != key:
            neighbors += [key]  # Convert x into a list containing x

        if isinstance(config['T0_cells'], int):
            config['T0_cells'] = [config['T0_cells']]  # Convert x into a list containing x
        if isinstance(config['T1_cells'], int):
            config['T1_cells'] = [config['T1_cells']]  # Convert x into a list containing x

        #cc = values_in_order[ids_in_order.index(key)]
        #fontColor = "white"
        #if cc > .33:
        #    fontColor = "black"
        #color = rgb_to_hex([cc,min(1.0,cc*3),min(1.0,cc*2)])
                           
        vv = states_in_order[ids_in_order.index(key)]
        fontColor = 'white'
        if vv == 0:
            color = rgb_to_hex([.2,.2,.2])
        if vv == -1:
            color = rgb_to_hex([1,1,0])
            fontColor = 'black'
        if vv > 0:
            color = rgb_to_hex([1,0,0])
        
            
        if key in config['T0_cells']:
            graph_lines.append(str(key) + ' [style = "filled",fillcolor = "' + color + '",fontcolor  = "'+fontColor+'", color="blue", penwidth=20]')
        elif key in config['T1_cells']:
            graph_lines.append(str(key) + ' [style = "filled",fillcolor = "' + color + '",fontcolor  = "'+fontColor+'", color="orange", penwidth=20]')
        else:   
            graph_lines.append(str(key) + ' [style = "filled",fillcolor = "' + color + '",fontcolor  = "'+fontColor+'"]')
        for i, neighbor in enumerate(neighbors):
            # Add edge from key to neighbor
            # load this cells predictions!

            this_fragMap=load_frag_map(inPath+'vis_output/cellPredictions/cellPrediction_'+str(key)+'.py')

            #if step == None:
            #    this_fragMap=load_frag_map(path+'vis_output/cellPredictions/cellPrediction_'+str(key)+'.py')
            #else:
            #    this_fragMap=load_frag_map(path+'vis_output/cellPredictions/cellPrediction_'+str(key)+'_step_'+str(step)+'.py')

            
            #graph_lines.append(f"    {key} -> {neighbor};")
            penWidth = this_fragMap[0][i] * connectionWeighting
            penWidth = this_fragMap[0][i] * connectionWeighting
            if penWidth < .2:
                penWidth = 0
            else:
                penWidth = 3 + 2 * penWidth
            if penWidth > 100:
                penWidth = 100


            #graph_lines.append(f'    {neighbor} -> {key} [penwidth={penWidth}, color="black", arrowsize = '+str(penWidth*arrowSize)+'];')
            graph_lines.append(f'    {neighbor} -> {key} [penwidth={penWidth}, color="black", arrowsize = '+str(.5+.125*penWidth*arrowSize)+'];')

            #print(graph_lines[-1])
            # Add edge from neighbor to key
            #graph_lines.append(f"    {neighbor} -> {key};")
        ## # add self
        ## penWidth = min(this_fragMap[0][len(neighbors)] * 5,100)
        ## if penWidth < .001:
        ##     penWidth = 0
        ## graph_lines.append(f'    {key} -> {key} [penwidth={penWidth}, color="black"];')

    count += 1
    graph_lines.append("}")  # End of the Graphviz DOT format
    return "\n".join(graph_lines)




###############################################################################################################
## FRAGMAP AND CONNECTOME TOOLS
###############################################################################################################

## first load the focal and neighbor IDs
def load_neighbor_ids(file_path):
    neighbor_dict = {}
    ids_in_order = []
    with open(file_path, 'r') as file:
        for line in file:
            key, value = line.strip().split(':')
            ids_in_order.append(int(key))
            neighbors = list(map(int, value.split(',')))
            neighbor_dict[int(key)] = neighbors
    return neighbor_dict,ids_in_order

def create_graphviz_text(neighbor_dict):
    graph_lines = ["digraph G {"]  # Start of the Graphviz DOT format
    for key, neighbors in neighbor_dict.items():
        for neighbor in neighbors:
            # Add edge from key to neighbor
            graph_lines.append(f"    {key} -> {neighbor};")
            # Add edge from neighbor to key
            #graph_lines.append(f"    {neighbor} -> {key};")
    graph_lines.append("}")  # End of the Graphviz DOT format
    return "\n".join(graph_lines)

def load_frag_map(file_name):
    frag_map_str = ''
    start_reading = False
    bracket_count = 0  # To count the opening and closing brackets
    with open(file_name, 'r') as file:
        for line in file:
            if 'fragMap' in line:
                start_reading = True
            if start_reading:
                # Count brackets to ensure completeness
                bracket_count += line.count('[')
                bracket_count -= line.count(']')
                frag_map_str += line
                # Only stop reading if all brackets are closed
                if bracket_count == 0:
                    break
    
    # Debug output to check the extracted string
    #print("Fragment Map String Extracted:", frag_map_str)

    # Locate the assignment to fragMap
    frag_map_start = frag_map_str.find('[')
    frag_map_end = frag_map_str.rfind(']') + 1
    frag_map_str = frag_map_str[frag_map_start:frag_map_end]

    # Safely evaluate the Python literal
    try:
        frag_map = ast.literal_eval(frag_map_str)
        return frag_map
    except SyntaxError as e:
        print("Error evaluating the fragment map:", e)
        print("Fragment Map String:", frag_map_str)
        return None  # Return None or handle the error appropriately

inPath = './'                                 #where to find files generated by MABE
tempPath = 'vis_output/'                      # a work space to store temp text and image files - well be cleaned up

outPath = 'images/'+str(rep)+'/'    # the directory related to each rep where final images are saved
path = os.path.join("images", str(rep))
# Create the directory if it doesn't exist
os.makedirs(path, exist_ok=True)

config_filename = 'config.csv'
config = load_csv_with_lists(inPath+config_filename)
config = convert_data(config)
print("loaded config:")
for key in config:
    if len(str(config[key])) > 30:
        print('   '+key+':',str(config[key])[:29],'...') # we don't need to see the entier genome!
    else:
        print('   '+key+':',config[key])

# Load some data from the file
neighbor_ids,ids_in_order = load_neighbor_ids(inPath+'vis_output/neighbor_ids.txt')

fragMap=load_frag_map(inPath+'vis_output/rawR_FragmentationMatrix.py')
f_count = len(fragMap) - 1 # this is the number of features in the fragmap map (in this case t-0 -> t-4)

# this takes the count of size 1 and 2 partions and returns the number of size 1 partitions (i.e., number of orig elements) -1 to make zero indexed
p1_count = int(round((1 + math.sqrt(1 + 8 * (len(fragMap[0])-2))) / 2)) - 1

print('partition size 1 count: ' + str(p1_count))
print('number of features: ' + str(f_count))

# now we will load all the step frag maps, we may not plot them, but loading them is fast and may be needed later
# outter list is the brain time steps while the inner lists are the features
step_FragMaps = []
for i in range(12):
    step_FragMaps.append(load_frag_map(inPath+'vis_output/rawR_step_'+str(i)+'_FragmentationMatrix.py'))


if renderFragMaps:
    ################################################################################
    #### FULL rawR fragmentaion matrix
    ################################################################################
    
    img = [l[:p1_count] for l in fragMap[:-1]]
    plt.figure(figsize=(10, 10))  # Width = 10 inches, Height = 5 inches
    plt.imshow(img,cmap='gray',aspect=5,interpolation='nearest',vmax=1.0)
    plt.yticks(range(len(features)),features)
    plt.xticks(range(len(ids_in_order)),[str(id) for id in ids_in_order],rotation = 90)
    plt.title("rawR for all paritions of size 1")
    
    plt.savefig(outPath+'/rawR_FULL.png')
    
    #plt.show()
    plt.close()

    ################################################################################
    #### FULL rawR fragmentaion matrix MOVE
    ################################################################################
    
    fragMap_MOVE=load_frag_map(inPath+'vis_output/rawR_FragmentationMatrix_MOVE.py')
    img = [l[:p1_count] for l in fragMap_MOVE[:-1]]
    plt.figure(figsize=(10, 10))  # Width = 10 inches, Height = 5 inches
    plt.imshow(img,cmap='gray',aspect=5,interpolation='nearest',vmax=1.0)
    plt.yticks(range(len(features)),features)
    plt.xticks(range(len(ids_in_order)),[str(id) for id in ids_in_order],rotation = 90)
    plt.title("** MOVE **  rawR for all paritions of size 1")
    
    plt.savefig(outPath+'/rawR_FULL_MOVE.png')
    
    #plt.show()
    plt.close()
    
    
    ################################################################################
    #### rawR STEPS (decompsed by brainStep) size 1 partitions fragmentaion matrix
    ################################################################################
    for j in range(f_count):      # nback offset
        img = []
        yTick_names = []
        for i in range(12):
            img.append([l for l in step_FragMaps[i][j][:p1_count]])
            print('step'+str(i)+': '+str(max(img[-1][:-3]))+ '  sum: ' + str(sum(img[-1][:-3])))
            yTick_names.append(features[j]+'\nstep '+str(i))
        
        plt.figure(figsize=(10, 10))  # Width = 10 inches, Height = 5 inches
        plt.imshow(img,cmap='gray',aspect=5,interpolation='nearest',vmax=1.0)
        plt.yticks(range(len(yTick_names)),yTick_names)
        plt.xticks(range(len(ids_in_order)),[str(id) for id in ids_in_order],rotation = 90)
        plt.title("rawR for all paritions of size 1 decomposed by brain update step")
        
        plt.savefig(outPath+'/rawR_parts_1_back'+str(j)+'.png')
        
        #plt.show()
        plt.close()

    ################################################################################
    #### rawR STEPS    MOVE     (decompsed by brainStep) size 1 partitions fragmentaion matrix
    ################################################################################
    step_FragMaps_MOVE = []
    for i in range(12):
        step_FragMaps_MOVE.append(load_frag_map(inPath+'vis_output/rawR_step_'+str(i)+'_FragmentationMatrix_MOVE.py'))
    for j in range(f_count):      # nback offset
        img = []
        yTick_names = []
        for i in range(12):
            img.append([l for l in step_FragMaps_MOVE[i][j][:p1_count]])
            print('step'+str(i)+': '+str(max(img[-1][:-3]))+ '  sum: ' + str(sum(img[-1][:-3])))
            yTick_names.append(features[j]+'\nstep '+str(i))
        
        plt.figure(figsize=(10, 10))  # Width = 10 inches, Height = 5 inches
        plt.imshow(img,cmap='gray',aspect=5,interpolation='nearest',vmax=1.0)
        plt.yticks(range(len(yTick_names)),yTick_names)
        plt.xticks(range(len(ids_in_order)),[str(id) for id in ids_in_order],rotation = 90)
        plt.title("** MOVE **  rawR for all paritions of size 1 decomposed by brain update step")
        
        plt.savefig(outPath+'/rawR_parts_1_back'+str(j)+'_MOVE.png')
        
        #plt.show()
        plt.close()

else:
    print("rendering fragMaps off")

##################################################################################################
# render state prediction per node (corrilated with world state/s)
##################################################################################################
if renderFullDataFlow:
    print("generating full data flow")
    # Generate the Graphviz DOT format text
    create_graphviz_rich(neighbor_ids,fragMap,ids_in_order,config,features)
else:
    print("full data flow off")


##################################################################################################
# render information flow graphs for each brain step (corrilated with world state/s)
##################################################################################################

if renderStepsDataFlow:
    print("generating steps data flow")
    create_graphviz_rich_steps(neighbor_ids,step_FragMaps,ids_in_order,config,features,timeColorConnections=True)
    create_graphviz_rich_steps(neighbor_ids,step_FragMaps,ids_in_order,config,features,timeColorConnections=False)
else:
    print("steps data flow off")

##################################################################################################
# convert text file from mabe to a state to state graph
##################################################################################################
if renderState2State:
    print('rendering state to state plot')
    convert_dot_txt_to_png(inPath+'vis_output/state_2_state.txt',outPath+'state_2_state.png')
else:
    print("render state2state off")






##################################################################################################
# create activation patterns and activation map
##################################################################################################
print("generating activation patterns")

def read_and_convert(filename):
    list_of_lists = []
    with open(filename, 'r') as file:
        for line in file:
            inner_list = []
            i = 0
            while i < len(line):
                if line[i] == '-':
                    if i+1 < len(line) and line[i+1] == '1':
                        inner_list.append(1)
                        i += 2  # Skip the '1' after '-'
                        continue
                elif line[i] == '0':
                    inner_list.append(0)
                elif line[i] == '1':
                    inner_list.append(.05)
                i += 1
            list_of_lists.append(inner_list)
    return list_of_lists

# find the line in the file with "active"; this gets the index in brain states for node "active"
def find_line_with_active(filename, active):
    active_prefix = str(active) + ":"  # Construct the prefix to match
    with open(filename, 'r') as file:
        for line_number, line in enumerate(file, start=1):
            if line.startswith(active_prefix):
                return line_number
    return -1  # Return -1 if no match is found

# see if two nodes have the "same" pattern
def compare_sublists(sublist1, sublist2, max_shift=24, skip=240):
    # Trim the first 50 elements
    sublist1 = sublist1[skip:]
    sublist2 = sublist2[skip:]

    length1 = len(sublist1)
    length2 = len(sublist2)

    best_ratio = 0

    # Compare sublists with shifts from -max_shift to max_shift
    for shift in range(-max_shift, max_shift + 1):
        score = 0
        # Calculate the indices to compare based on shift
        if shift >= 0:
            # Shifting sublist2 to the right
            start1 = 0
            start2 = shift
            length_to_compare = min(length1, length2 - shift)
        else:
            # Shifting sublist2 to the left
            start1 = -shift
            start2 = 0
            length_to_compare = min(length1 - start1, length2)

        # Calculate score for this alignment
        for i in range(length_to_compare):
            if sublist1[start1 + i] == sublist2[start2 + i]:
                score += 1

        # Update best ratio if this one is better
        if length_to_compare > 0:
            current_ratio = score / length_to_compare
            best_ratio = max(best_ratio, current_ratio)

    return best_ratio

#Load brain states
filename = "vis_output/allBrainStates.txt"
brain_states = read_and_convert(filename)

def load_neighbor_ids(file_path):
    neighbor_dict = {}
    ids_in_order = []
    with open(file_path, 'r') as file:
        for line in file:
            key, value = line.strip().split(':')
            ids_in_order.append(int(key))
            neighbors = list(map(int, value.split(',')))
            neighbor_dict[int(key)] = neighbors
    return neighbor_dict,ids_in_order

path = './'
neighbor_ids,nodesList = load_neighbor_ids(path+'vis_output/neighbor_ids.txt')

# make an list of list where each brains pattern is broken out
allActivations = [[] for _ in nodesList]
for brain_state in brain_states[:12*180]:
    for NINDEX in range(len(nodesList)):
        allActivations[NINDEX].append(brain_state[NINDEX])

# for each NID in nodesList check agaist every type in types Dict
# first NID is simply inserted
# we will then have a Dict where each key contains a list of the nodes with the same signal pattern!

typesIndex = 0
typesDict = {}
typesDict[typesIndex] = [nodesList[0]]
typesIndex += 1

for NINDEX in range(1, len(nodesList)):
    NewType = True
    bestMatch = 0.0
    for typeKey in typesDict.keys():
        score = compare_sublists(allActivations[NINDEX], allActivations[nodesList.index(typesDict[typeKey][0])], max_shift=24, skip=12*20)
        bestMatch = max(bestMatch,score)
        #print(f"Best alignment score (ratio)  node: node{nodesList[NINDEX]} / type_{typeKey} = {score}")
        
        if score == 1.0:
            #print('    node:',nodesList[NINDEX],' is type',typeKey)
            typesDict[typeKey].append(nodesList[NINDEX])
            NewType = False
            break  # This breaks out of the inner for-loop

    if NewType:
        typesDict[typesIndex] = [nodesList[NINDEX]]
        #print('NODE:',nodesList[NINDEX],'    NEWTYPE:',typesIndex,'   with best match:',bestMatch)
        typesIndex += 1

colorList = [
    (255, 0, 0),       # Red
    (0, 255, 0),       # Lime
    (60, 60, 255),       # Blue (original, will be adjusted to a lighter shade below)
    (255, 255, 0),     # Yellow
    (0, 255, 255),     # Cyan
    (255, 0, 255),     # Magenta
    (192, 192, 192),   # Silver
    (188, 60, 60),       # Maroon
    (168, 168, 0),     # Olive
    (60, 188, 60),       # Green
    (178, 0, 178),     # Purple
    (0, 178, 178),     # Teal
    (135, 206, 235),   # Light Sky Blue (replacing the original Blue)
    (255, 165, 0),     # Orange
    (255, 20, 147),    # Deep Pink
    (165, 40, 180),      # Indigo
    (255, 192, 203),   # Pink
    (218, 165, 32),    # Goldenrod
    (106, 90, 205),    # Slate Blue
    (40, 120, 40)        # Dark Green
] * 1000

for i in range(len(colorList)):
    colorList[i] = (colorList[i][0]/255,colorList[i][1]/255,colorList[i][2]/255)

i = 0

jointFig = []
nodesList = [typesDict[key][0] for key in typesDict.keys()]

# Create the figure object to hold the subplots
fig = plt.figure(figsize=(len(nodesList)*1, 6))  # Adjust the figure size as needed

for active,LABEL in zip(nodesList,[str(n) for n in nodesList]):
    #print(i,active)
    LABEL = str(i)
    activeIndex = find_line_with_active("vis_output/neighbor_ids.txt", active) - 1
    #print(active,activeIndex)
    ax = fig.add_subplot(1, len(nodesList), i+1) # i+1 indicates the current subplot
    
    img = []
    line = []
    for brain_state in brain_states:
        if len(line) == 12:
            img.append(line)
            line = []
        line.append(brain_state[activeIndex])
        
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", [(0,0,0),colorList[i]])

    ax.imshow(img[0:50],cmap=custom_cmap,aspect=1,interpolation="NONE")
    for x in np.arange(.5,12,1):
        plt.axvline(x=x, color=(.6,.6,.6), linestyle='-',linewidth=1)  # Customize color and linestyle as needed
    #ax.set_xticks(range(0,12),[str(x) for x in range(0,12)],rotation=90)
    ax.set_xlabel('step')
    if i == 0:
        ax.set_ylabel('world update')
        ax.set_yticks([0,9,19,29,39,49])
        ax.set_yticklabels(['1','10','20','30','40','50'])
    if i > 0:  # Condition to hide the y-axis for all subplots except the first
        ax.set_yticks([])  # This hides the y-ticks themselves
        ax.set_yticklabels([])  # This hides the y-tick labels

    ax.set_title(str(i),fontsize=15)
    ax.set_xticks([0,5,11])
    ax.set_xticklabels(['1','6','12'])
    i+=1
plt.tight_layout()
plt.savefig(outPath+'/activationPatterns.png')
#plt.show()
plt.close()



#### MAKE THE GRAPH / ACTIVATION MAP
print("generating activation map")

def find_key_by_value(my_dict, value):
    for key, values in my_dict.items():
        if value in values:
            return key
    return None  # Return None if the value is not found in any list
def rgb_to_hex(rgb):
    """
    Convert an RGB list [R, G, B] where each value is between 0 and 1
    into a hexadecimal string '#XXXXXX'.
    """
    return '#{:02X}{:02X}{:02X}'.format(
        int(rgb[0] * 255), 
        int(rgb[1] * 255), 
        int(rgb[2] * 255)
    )

neighbor_ids,nodesList = load_neighbor_ids(path+'vis_output/neighbor_ids.txt')

graph_lines = ["digraph G {"]  # Start of the Graphviz DOT format
graph_lines.append('graph [ranksep=0.125];')
graph_lines.append('node [shape=ellipse, fixedsize=true, width=1, height=.7,fontname="Times-Bold"];')

allLinks = []
arrowSize = .4
penWidth = 5

for NID in nodesList:
    key = find_key_by_value(typesDict,NID)
    #key = nodesList.index(typesDict[key][0])
    color = colorList[key]
    graph_lines.append(str(NID)+' [label="'+str(key)+'",style="filled",fillcolor="'+rgb_to_hex(color)+'", color="black", penwidth=3,fontsize=30]')
    for neighbor in neighbor_ids[NID]:
        if str(neighbor)+'->'+str(NID) in allLinks:
            graph_lines.append('    '+str(NID)+' -> '+str(neighbor)+' [penwidth=0, arrowsize = 0]')
        else:
            graph_lines.append('    '+str(NID)+' -> '+str(neighbor)+' [dir=both,penwidth='+str(penWidth)+', color="black", arrowsize = 0]')
            allLinks.append(str(NID)+'->'+str(neighbor));
graph_lines.append('}')
graphviz_text = "\n".join(graph_lines)

tempPath = 'vis_output/'
tempFile = tempPath+'temp.txt'               # name of temp file with graphviz code
with open(tempFile, 'w') as file:
    file.write(graphviz_text)
convert_dot_txt_to_png(tempFile,outPath+'/activationMap.png')







##################################################################################################
# render brain states on flow graph to make brain playback frames
##################################################################################################

def parse_line(line):
    # Initialize an empty list to collect parsed integers
    parsed = []
    i = 0
    while i < len(line):
        if line[i] == '-':
            if i + 1 < len(line) and line[i+1] == '1':
                # Append -1 to the list and skip the next character
                parsed.append(-1)
                i += 2  # move past "-1"
                continue
        elif line[i] in '01':
            # Append 0 or 1 as integers
            parsed.append(int(line[i]))
        i += 1
    return parsed

def read_brain_states(file_path):
    brain_states = []
    with open(file_path, 'r') as file:
        for line in file:
            # Parse each line and append the result to the brain_states list
            parsed_line = parse_line(line.strip())
            brain_states.append(parsed_line)
    return brain_states

file_path = inPath+'vis_output/allBrainStates.txt'
brain_state_list = read_brain_states(file_path)

if(renderUpdates > 0):
    print('rendering updates ',end='',flush=True)

    # Save or display the modified image
    if not os.path.exists(outPath+'updates'):
        os.makedirs(outPath+'updates')

    for u in range(renderUpdates):
        print('.',end='',flush=True)
        if u%12 == 0:
            bg_color = '#FFFFEE'
        else:
            bg_color = 'white'
        # Generate the Graphviz DOT format text
        graphviz_text = create_graphviz_brainState(neighbor_ids,[l[:p1_count] for l in fragMap[:-1]][0],brain_state_list[u],ids_in_order,config)   #,step=u%12)
        #    file.write(graphviz_text)
        tempFile = tempPath+'temp.txt'               # name of temp file with graphviz code
        outFile = f"{outPath}updates/update_{u:04}.png" # name of file with rendered graph

        with open(tempFile, 'w') as file:
            file.write(graphviz_text)
        convert_dot_txt_to_png(tempFile,outFile)
        # Load the image    
        img = Image.open(outFile)

        # Specify the height of the additional space for text
        additional_height = 400  # Adjust as needed for the amount of text
        
        # Create a new image with extra space at the bottom
        new_img = Image.new('RGB', (img.width, img.height + additional_height), (255, 255, 255))
        new_img.paste(img, (0, 0))

        # Prepare to draw on the image
        draw = ImageDraw.Draw(new_img)
        
        # Optional: Load a font (or use a default one)
        font = ImageFont.truetype("arial.ttf", size=50)  # You can specify your own font
        #font = ImageFont.load_default()
        
        # Define the text and its position
        if u%12 == 0:
            text = "Update: " +str(u//12) + "  Frame: " + str(u) + "("+str(u%12)+")   **World Update**\n"
        else:
            text = "Update: " +str(u//12) + "  Frame: " + str(u) + "("+str(u%12)+")\n"
        text += "sensor input: " + str(config['T0_cells'])+"\n\n"
        text += "yellow.... charge\n"
        text += "red....... reset\n"
        position = (10, new_img.height - 300)  # Adjust position as needed
        draw.text(position, text, font=font, size=50, fill="black")

        text = "\n"
        text += "sensor wall bump: " + str(config['T0_cells'][0])+"\n"
        text += "sensor wall front: " + str(config['T0_cells'][1])+"\n"
        text += "sensor direction A: " + str(config['T0_cells'][2])+"\n"
        text += "sensor direction B: " + str(config['T0_cells'][3])+"\n"

        position = (new_img.width * .25, new_img.height - 300)  # Adjust position as needed
        draw.text(position, text, font=font, size=50, fill="black")

        new_img.save(outFile)
        new_img.close()
    print()
    print('done')
else:
    print("rendering updates off")
    
    