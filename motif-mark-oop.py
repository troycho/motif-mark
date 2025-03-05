#!/usr/bin/env python3

import argparse
import re
import random
import cairo


# Global variables for Pycairo 
LEFT_PADDING = 20 # padding in left margin
LEGEND_PADDING = 40  # Space between last gene group and legend
LEGEND_ITEM_SPACING = 40  # Space between each legend item
LEGEND_RECT_WIDTH = 30  # Width of rectangles in legend
LEGEND_RECT_HEIGHT = 10  # Height of rectangles in legend
TEXT_OFFSET = 3  # Space between rectangle and text

class Motif:

    MOTIF_HEIGHT = 10  # Height of motif rectangle
    TRANSPARENCY = 0.7  # Set transparency level (0 = fully transparent, 1 = fully opaque)

    def __init__(self, start, stop, color):
        self.start = start
        self.stop = stop
        self.color = color  # (R, G, B) tuple

    def draw(self, ctx, y_position):
        ctx.set_antialias(cairo.ANTIALIAS_NONE)  
        
        # Use RGBA (adds transparency)
        ctx.set_source_rgba(self.color[0], self.color[1], self.color[2], self.TRANSPARENCY)  

        # Draw the motif rectangle
        ctx.rectangle(self.start + LEFT_PADDING, y_position - self.MOTIF_HEIGHT / 2, 
                      self.stop - self.start, self.MOTIF_HEIGHT)
        ctx.fill()

        ctx.set_antialias(cairo.ANTIALIAS_DEFAULT)  # Reset anti-aliasing
    
    def __repr__(self):
        return f"Color: {self.color}, Position: ({self.start}, {self.stop})"

class Exon:

    EXON_HEIGHT = 10
    EXON_COLOR = (0, 0, 0)  # Black

    def __init__(self, start, stop):
        self.start = start
        self.stop = stop
    
    def draw(self, ctx, y_position):
        ctx.set_antialias(cairo.ANTIALIAS_NONE)  # Disable anti-aliasing
        ctx.set_source_rgb(*self.EXON_COLOR)  # Set fill color
        ctx.rectangle(self.start + LEFT_PADDING, y_position - self.EXON_HEIGHT / 2, 
                      self.stop - self.start, self.EXON_HEIGHT)
        ctx.fill()
        ctx.set_antialias(cairo.ANTIALIAS_DEFAULT)  # Reset anti-aliasing to default


    def __repr__(self):
        return f"Exstart: {self.start}, Exstop: {self.stop}"

class Gene:
    def __init__(self, start, stop, name):
        self.start = start
        self.stop = stop
        self.name = name

    def draw(self, ctx, y_position):
        ctx.set_antialias(cairo.ANTIALIAS_NONE)
        ctx.set_source_rgb(0, 0, 0)  # Black line
        ctx.set_line_width(2)
        ctx.move_to(self.start + LEFT_PADDING, y_position)
        ctx.line_to(self.stop + LEFT_PADDING, y_position)
        ctx.stroke()
        ctx.set_antialias(cairo.ANTIALIAS_DEFAULT)

        # Draw gene name above the line
        ctx.set_source_rgb(0, 0, 0)
        ctx.select_font_face("Arial") #cairo.FONT_SLANT_NORMAL
        ctx.set_font_size(18)
        ctx.move_to(self.start + LEFT_PADDING, y_position - 20)
        ctx.show_text(self.name)
    
    def __repr__(self):
        return f"Gene name: {self.name}"

class GeneGroup:
    # global class variable that should store all instantiated objects
    GENEGROUP_LIST = []

    def __init__(self, exon: Exon, motifs: list[Motif], gene: Gene):
        self.exon = exon
        self.motifs = motifs
        self.gene = gene
        GeneGroup.GENEGROUP_LIST.append(self)

    def draw(self, ctx, y_position):
        self.gene.draw(ctx, y_position)
        self.exon.draw(ctx, y_position)
        for motif in self.motifs:
            motif.draw(ctx, y_position)

    def __repr__(self):
        return f"Genegroup = Exon: {self.exon}\nMotifs: {self.motifs}\nGene: {self.gene}"

def get_args():
    parser = argparse.ArgumentParser(description="Visualize protein binding motifs in gene sequences")
    parser.add_argument("-f", "--fasta", help="Input fasta file")
    parser.add_argument("-m", "--motifs", help="Input motifs file")
    return parser.parse_args()

def oneline_fasta(file: str, out_file: str) -> None:
    '''Takes path to FASTA file and output file name as arguments. Converts each record in file to single-line FASTA record and writes to out_file. 
    Returns None.'''
    with open(file, "r") as fh, open(out_file, "w") as wf:
        first_time: bool=True
        for line in fh:
            if first_time:
                wf.write(line)
                first_time = False
            elif ">" in line:
                wf.write(f'\n{line}')
            else:
                wf.write(line.strip("\n"))

def build_regex(motif: str) -> str:
    '''Returns a regex pattern that matches all possible combinations of
    the input motif sequence (which needs to be all lowercase).'''
    mapping = {'y': '[ctu]',
               't': '[ut]',
               'u': '[ut]',
               'r': '[ag]',
               'm': '[ac]',
               'k': '[gt]',
               's': '[gc]',
               'w': '[at]',
               'h': '[act]',
               'b': '[cgt]',
               'v': '[acg]',
               'd': '[agt]',
               'n': '[atcgu]'}

    # 'y' can be any pyrimidine base; keep all other bases the same
    pattern = ''.join(mapping.get(base, base) for base in motif)
    return pattern

def find_motif(seq: str, motif: str) -> list:
    '''Takes gene sequence and motif sequence and returns a list of tuples containing start and end indices of
    all sites in the sequence where motif was found.'''
    seq = seq.lower()
    motif = motif.lower()
    pattern = build_regex(motif) #get generalized regex pattern of motif
    matches: list = []
    for i in range(len(seq)):
        match = re.match(pattern, seq[i:])
        if match:
            matches.append((i + match.start(), i + match.end() - 1))
    return matches


if __name__ == "__main__":
    args = get_args()

    # UNCOMMENT THE LINE BELOW BEFORE FINAL SUBMISSION
    oneline_fasta(args.fasta, f"oneline_{args.fasta}")

    # 5 RGB colors
    color_set = {(0.055, 0.722, 0.733),
                 (0.471, 0.369, 0.941),
                 (0.863, 0.149, 0.498),
                 (0.996, 0.380, 0),
                 (1, 0.690, 0)}
    
    # Assign a different color code to each motif
    motif_dict: dict = {}

    with open(args.motifs, "r") as fh:
        for line in fh:
            line = line.strip()
            motif_dict[line] = random.choice(list(color_set))
            color_set.remove(motif_dict[line])

    with open(f"oneline_{args.fasta}", "r") as fh:
        record: list = [] # to hold current record
        for i, line in enumerate(fh):
            line = line.strip()
            record.append(line)
            if i % 2 != 0:
                motif_objects = []
                for motif in motif_dict.keys():
                    positions = find_motif(record[1], motif) #list of tuples of indices where motif occurs in sequence
                    for pos in positions: #create unique Motif object for each occurrence
                        mo = Motif(pos[0], pos[1], motif_dict[motif])
                        motif_objects.append(mo)
                exon = None # will be assigned to Exon object
                try:
                    exon_match = re.search("[ACTGU]+", record[1])
                    exon = Exon(exon_match.start(), exon_match.end()-1)
                except Exception as e:
                    print("There's no exon in this gene, that can't be right.")
                gene = Gene(0, len(record[1])-1, record[0])
                gene_group = GeneGroup(exon, motif_objects, gene)
                record.clear()

    # Pycairo figure
    WIDTH = 1000
    HEIGHT = 100 + len(GeneGroup.GENEGROUP_LIST) * 70  # Dynamic height
    
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)

    # Background color
    ctx.set_source_rgb(1, 1, 1)
    ctx.paint()

    # Padding and positioning
    x_padding = 50
    y_padding = 50
    y_pos = 0

    # Scale factor (adjust this to fit the longest sequence in width)
    max_gene_length = max(gene_group.gene.stop for gene_group in GeneGroup.GENEGROUP_LIST)
    scale_factor = (WIDTH - 2 * x_padding) / max_gene_length

    # Draw each GeneGroup
    for i, gene_group in enumerate(GeneGroup.GENEGROUP_LIST):
        y_pos = y_padding + i * 80
        gene_group.gene.start *= scale_factor
        gene_group.gene.stop *= scale_factor
        gene_group.exon.start *= scale_factor
        gene_group.exon.stop *= scale_factor
        for motif in gene_group.motifs:
            motif.start *= scale_factor
            motif.stop *= scale_factor

        gene_group.draw(ctx, y_pos)

    # Draw legend
    x_start = 5 
    y_start = y_pos + 50 

    ITEM_SPACING = 130  # Space between legend elements

    # Gene legend
    ctx.set_source_rgb(0, 0, 0)
    ctx.move_to(x_start, y_start)
    ctx.line_to(x_start + LEGEND_RECT_WIDTH, y_start)
    ctx.stroke()
    ctx.move_to(x_start + LEGEND_RECT_WIDTH + TEXT_OFFSET, y_start + 3)
    ctx.show_text("Gene")

    x_start += 100  # Move right for next item

    # Exon legend
    ctx.set_source_rgb(0, 0, 0)
    ctx.rectangle(x_start, y_start - LEGEND_RECT_HEIGHT / 2, LEGEND_RECT_WIDTH, LEGEND_RECT_HEIGHT)
    ctx.fill()
    ctx.move_to(x_start + LEGEND_RECT_WIDTH + TEXT_OFFSET, y_start + 3)
    ctx.show_text("Exon")

    x_start += ITEM_SPACING  # Move right for next item

    # Draw Motif Legends (each motif horizontally aligned)
    for motif, color in motif_dict.items():
        ctx.set_source_rgb(*color)  # Use motif color
        ctx.rectangle(x_start, y_start - LEGEND_RECT_HEIGHT / 2, LEGEND_RECT_WIDTH, LEGEND_RECT_HEIGHT)
        ctx.fill()
        ctx.move_to(x_start + LEGEND_RECT_WIDTH + TEXT_OFFSET, y_start + 3)
        ctx.show_text(motif)

        x_start += ITEM_SPACING  # Move right for the next item

    # Save the figure
    prefix = args.fasta.split(".")
    surface.write_to_png(f"{prefix[0]}.png")
