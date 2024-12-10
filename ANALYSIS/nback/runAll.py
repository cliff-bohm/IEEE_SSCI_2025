import matplotlib.pyplot as plt
import math
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
renderUpdates = 23*12#27*12 # 240 # 2400 per lifetime!
print("renderUpdates = ",renderUpdates)

features = ['t-0','t-1','t-2','t-3','t-4']
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

        text = "\nt-1 output: " + str(config['T1_cells'][0])+"\n"
        text += "t-2 output: " + str(config['T1_cells'][1])+"\n"
        text += "t-3 output: " + str(config['T1_cells'][2])+"\n"
        text += "t-4 output: " + str(config['T1_cells'][3])+"\n"
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

        text = "\nt-1 output: " + str(config['T1_cells'][0])+"\n"
        text += "t-2 output: " + str(config['T1_cells'][1])+"\n"
        text += "t-3 output: " + str(config['T1_cells'][2])+"\n"
        text += "t-4 output: " + str(config['T1_cells'][3])+"\n"
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
    
    allLinks = []
    
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
            graph_lines.append(str(key) + ' [style = "filled",fillcolor = "' + color + '",fontcolor  = "'+fontColor+'", color="black", penwidth=3]')
        for i, neighbor in enumerate(neighbors[:-1]):
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

            penWidth = 8 # hardcoded for now!
            if str(key)+'->'+str(neighbor) in allLinks:
                graph_lines.append(f'    {neighbor} -> {key} [penwidth=0, arrowsize = 0];')
            else:
                #graph_lines.append(f'    {neighbor} -> {key} [dir=both,penwidth={penWidth}, color="black", arrowsize = '+str(.5+.125*penWidth*arrowSize)+'];')
                graph_lines.append(f'    {neighbor} -> {key} [dir=both,penwidth={penWidth}, color="black", arrowsize = 0];')
                allLinks.append(str(neighbor)+'->'+str(key));
            
            graph_lines.append(f'    {key} -> {key} [penwidth=0, arrowsize = 0];') # include self link for spacing

            #graph_lines.append(f'    {neighbor} -> {key} [penwidth={penWidth}, color="black", arrowsize = '+str(.5+.125*penWidth*arrowSize)+'];')

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

    for u in range(36,renderUpdates):
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

        text = "\nt-1 output: " + str(config['T1_cells'][0])+"\n"
        text += "t-2 output: " + str(config['T1_cells'][1])+"\n"
        text += "t-3 output: " + str(config['T1_cells'][2])+"\n"
        text += "t-4 output: " + str(config['T1_cells'][3])+"\n"
        position = (new_img.width * .25, new_img.height - 300)  # Adjust position as needed
        draw.text(position, text, font=font, size=50, fill="black")

        new_img.save(outFile)
        new_img.close()
    print()
    print('done')
else:
    print("rendering updates off")
