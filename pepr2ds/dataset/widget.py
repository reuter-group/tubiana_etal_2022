from .tagibs import *

import ipywidgets as widgets
import nglview as nv





class Widgets():
        def __init__(self, parent):
            self.parent = parent
            self.create_widgets()
            self.watch_widgets()
            self.set_layout()

            #instantiate features
            self.alignment_group = []

        def create_widgets(self):
            self.dropdownDomainsList = self.parent.dataset.domain.unique()  # drodown
            self.domainsWidget = widgets.Dropdown(description="Domain",
                                                  options=self.dropdownDomainsList,
                                                  value='PH')
            self.pdbWidget = widgets.Dropdown(description='PDB Ref',
                                              disabled=True)
            self.extendSSWidget = widgets.Checkbox(value=True,
                                                   description='Extend SS ?')
            self.alignmentWidget = widgets.Checkbox(value=True,
                                                    description='enrich with 2D prosite sequence ?')
            self.onlyCoilSSWidget = widgets.Checkbox(value=False,
                                                     description='Only Coil ?')
            self.cathClusterWidget = widgets.Dropdown(description="Cath Cluster",
                                                      options=[None, "S35", "S60", "S95", "S100"],
                                                      value="S100")
            self.showCoinsertableWidget = widgets.Checkbox(value=False,
                                                           description='Show co-insertables ?')

            self.showProtrusionWidget = widgets.Checkbox(value=False,
                                                         description='Show Protrusion ?')
            self.excludeStrand = widgets.Checkbox(value=False,
                                                  description='ExcludeStrand ?')

            self.unirefWidget = widgets.Dropdown(description="Uniref",
                                                 options=[None, "uniref50", "uniref90", "uniref100"],
                                                 value="uniref100")

            self.IncludeFromWidget = widgets.IntText(value=20, description='From:', disabled=True, )
            self.IncludeToWidget = widgets.IntText(value=26, description='To:', disabled=True, )
            self.IncludeSegments = [widgets.HBox(children=[self.IncludeFromWidget, self.IncludeToWidget])]

            self.ExcludeFromWidget = widgets.IntText(value=0, description='From:', disabled=True, )
            self.ExcludeToWidget = widgets.IntText(value=0, description='To:', disabled=True, )
            self.ExcludeSegments = [widgets.HBox(children=[self.ExcludeFromWidget, self.ExcludeToWidget])]

            self.entropyWidget = widgets.FloatSlider(value=4.22, min=0, max=4.22,
                                                     step=0.01, description='', )
            self.expositionWidget = widgets.IntSlider(value=100, min=0, max=100,
                                                      step=1, description='', )

            self.entropySignWidget = widgets.Dropdown(description="Entropy",
                                                      options=['<', '<=', '=', '>', '>='],
                                                      value='<=',
                                                      layout=widgets.Layout(width="20%"))

            self.expositionSignWidget = widgets.Dropdown(description="Exposition (%)",
                                                         options=['<', '<=', '=', '>', '>='],
                                                         value='<=',
                                                         layout=widgets.Layout(width="20%"))

            # Buttons
            self.createDatasetButton = widgets.Button(description="Create subdataset")
            self.residueAnalysisButton = widgets.Button(description="Amino acid Analysis", disabled=True)
            self.addIncludeSegmentButton = widgets.Button(description="Add new Segment", disabled=True)
            self.addExcludeSegmentButton = widgets.Button(description="Add new Segment", disabled=True)
            self.checkSelectionButton = widgets.Button(description="Create ConvexHull", disabled=True)
            self.quickRunPHButton = widgets.Button(description="Quick setup PH")
            self.quickRunC2Button = widgets.Button(description="Quick setup C2")

            # Output
            self.out = widgets.Output()
            self.displayMol = nv.NGLWidget()  # widgets.Text("select a PDB")

        def watch_widgets(self):
            # Action
            self.domainsWidget.observe(self.update_domainsWidget, names=['value'])
            self.pdbWidget.observe(self.update_pdbWidget, names=['value'])
            self.createDatasetButton.on_click(self.on_run_button_clicked)
            self.residueAnalysisButton.on_click(self.on_residueAnalysis_button_clicked)
            self.addIncludeSegmentButton.on_click(self.on_addIncludeSegment_button_clicked)
            self.addExcludeSegmentButton.on_click(self.on_addExcludeSegment_button_clicked)
            self.checkSelectionButton.on_click(self.on_createConvexhull_button_clicked)
            self.checkSelectionButton.on_click(self.update_selection)
            self.expositionWidget.observe(self.update_selection, names=["value"])
            self.entropyWidget.observe(self.update_selection, names=["value"])
            self.IncludeFromWidget.observe(self.update_selection, names=["value"])
            self.IncludeToWidget.observe(self.update_selection, names=["value"])
            self.ExcludeFromWidget.observe(self.update_selection, names=["value"])
            self.ExcludeToWidget.observe(self.update_selection, names=["value"])
            self.showCoinsertableWidget.observe(self.on_coinsertable_ticked, names=["value"])
            self.showProtrusionWidget.observe(self.on_protrusion_ticked, names=["value"])
            self.quickRunPHButton.on_click(self.quickRunPH)
            self.quickRunC2Button.on_click(self.quickRunC2)

        def set_layout(self):
            self.tab1HBox_layout = widgets.Layout(display='flex',
                                                  # flex_flow='column',
                                                  align_items='stretch',
                                                  border='solid',
                                                  width='90%')

            self.IncludeVBoxSegments = widgets.VBox(children=self.IncludeSegments)
            self.MetaIncludeVboxSegment = widgets.VBox(
                children=[self.IncludeVBoxSegments, self.addIncludeSegmentButton])

            self.ExcludeVBoxSegments = widgets.VBox(children=self.ExcludeSegments)
            self.ExcludeMetaVboxSegment = widgets.VBox(
                children=[self.ExcludeVBoxSegments, self.addExcludeSegmentButton])

            self.entropyBox = widgets.HBox([self.entropySignWidget, self.entropyWidget])
            self.expositionBox = widgets.HBox([self.expositionSignWidget, self.expositionWidget])

            self.selectionTab = widgets.Tab()
            self.selectionTab.set_title(0, 'Include in IBS')
            self.selectionTab.set_title(1, 'Exclude from IBS')
            self.selectionTab.children = [self.MetaIncludeVboxSegment, self.ExcludeMetaVboxSegment]
            # TAB1 BOX. Configuration
            self.tab1VBox = widgets.VBox(children=[widgets.HBox(children=[self.quickRunPHButton,
                                                                          self.quickRunC2Button]),
                                                   self.domainsWidget,
                                                   self.pdbWidget,
                                                   self.selectionTab,
                                                   self.entropyBox,
                                                   self.expositionBox,
                                                   self.showCoinsertableWidget,
                                                   self.showProtrusionWidget,
                                                   self.checkSelectionButton,
                                                   self.extendSSWidget,
                                                   self.alignmentWidget,
                                                   self.onlyCoilSSWidget,
                                                   self.excludeStrand,
                                                   self.cathClusterWidget,
                                                   self.unirefWidget,
                                                   self.createDatasetButton])
            self.tab1 = widgets.HBox(children=[self.tab1VBox,
                                               self.displayMol],
                                     layout=self.tab1HBox_layout)
            self.tab2 = widgets.VBox([self.residueAnalysisButton])

        def display(self):

            self.masterTab = widgets.Tab()
            self.masterTab.set_title(0, 'Create IBS dataset')
            self.masterTab.set_title(1, 'Analysis')

            self.masterTab.children = [self.tab1,
                                       self.tab2]
            display(self.masterTab)
            display(self.out)

            self.update_domainsWidget(self.domainsWidget.value)

        def quickRunPH(self, change):
            with self.out:
                self.domainsWidget.value = "PH"
                self.pdbWidget.value = "2da0A00"
                # Add sew segment
                newSegment1 = widgets.HBox([widgets.IntText(value=20, description='From:', disabled=False),
                                            widgets.IntText(value=26, description='to:', disabled=False)])
                newSegment1.children[0].observe(self.update_selection, names=["value"])
                newSegment1.children[1].observe(self.update_selection, names=["value"])

                newSegment2 = widgets.HBox([widgets.IntText(value=42, description='From:', disabled=False),
                                            widgets.IntText(value=50, description='to:', disabled=False)])
                newSegment2.children[0].observe(self.update_selection, names=["value"])
                newSegment2.children[1].observe(self.update_selection, names=["value"])

                self.IncludeSegments = [newSegment1, newSegment2]
                self.IncludeVBoxSegments.children = tuple(self.IncludeSegments)

                self.update_selection()
                self.unirefWidget.value = "uniref90"
                self.cathClusterWidget.value = "S95"

        def quickRunC2(self, change):
            with self.out:
                self.domainsWidget.value = "C2"
                self.pdbWidget.value = "1rsyA00"
                # Add sew segment
                newSegment1 = widgets.HBox([widgets.IntText(value=173, description='From:', disabled=False),
                                            widgets.IntText(value=177, description='to:', disabled=False)])
                newSegment1.children[0].observe(self.update_selection, names=["value"])
                newSegment1.children[1].observe(self.update_selection, names=["value"])

                newSegment2 = widgets.HBox([widgets.IntText(value=232, description='From:', disabled=False),
                                            widgets.IntText(value=238, description='to:', disabled=False)])
                newSegment2.children[0].observe(self.update_selection, names=["value"])
                newSegment2.children[1].observe(self.update_selection, names=["value"])

                self.IncludeSegments = [newSegment1, newSegment2]
                self.IncludeVBoxSegments.children = tuple(self.IncludeSegments)

                self.update_selection()
                self.unirefWidget.value = "uniref90"
                self.cathClusterWidget.value = "S95"

        def on_addIncludeSegment_button_clicked(self, change):
            """
            Add a new residue range.
            """
            with self.out:
                newSegment = widgets.HBox([widgets.IntText(value=0, description='From:', disabled=False),
                                           widgets.IntText(value=0, description='to:', disabled=False)])
                newSegment.children[0].observe(self.update_selection, names=["value"])

                newSegment.children[1].observe(self.update_selection, names=["value"])
                self.IncludeSegments += [newSegment]
                self.IncludeVBoxSegments.children = tuple(self.IncludeSegments)

        def on_addExcludeSegment_button_clicked(self, change):
            """
            Add a new residue range.
            """
            with self.out:
                newSegment = widgets.HBox([widgets.IntText(value=0, description='From:', disabled=False),
                                           widgets.IntText(value=0, description='to:', disabled=False)])
                newSegment.children[0].observe(self.update_selection, names=["value"])

                newSegment.children[1].observe(self.update_selection, names=["value"])
                self.ExcludeSegments += [newSegment]
                self.ExcludeVBoxSegments.children = tuple(self.ExcludeSegments)

        def update_domainsWidget(self, *args):
            domain = self.domainsWidget.value
            cathpdblist = self.parent.dataset.query("domain == @domain and data_type != 'prosite'").cathpdb.unique()
            cathpdblist = sorted(cathpdblist)
            self.pdbWidget.disabled = False
            self.pdbWidget.options = cathpdblist

        def update_pdbWidget(self, *args):
            with self.out:
                self.IncludeFromWidget.disabled = False
                self.IncludeToWidget.disabled = False
                self.ExcludeFromWidget.disabled = False
                self.ExcludeToWidget.disabled = False
                self.checkSelectionButton.disabled = False
                self.addIncludeSegmentButton.disabled = False
                self.addExcludeSegmentButton.disabled = False

                self.showCoinsertableWidget.value = False
                self.showProtrusionWidget.value = False
                domain = self.domainsWidget.value
                pdbname = self.pdbWidget.value

                if domain and pdbname:
                    structure = self.parent.CATHFOLDER + "domains/" + domain + "/zaligned/" + pdbname + ".pdb"
                    self.mdaMol = mda.Universe(structure)
                    # cog = self.mdaMol.atoms.center_of_geometry()
                    # self.mdaMol.atoms.translate(-cog)

                    #self.set_entropy_as_bfactor(pdbname)

                    self.displayMol = nv.show_mdanalysis(self.mdaMol)
                    self.displayMol.clear_representations()
                    self.displayMol.add_representation('cartoon', selection='protein', color='sstruc')
                    # self.displayMol.add_representation('cartoon', selection='protein',
                    #                                    color='bfactor')  # Representation 0 = Cartoon
                    self.displayMol.add_representation('hyperball', selection='@1',
                                                       color='red')  # Representation 1 = selected residues
                    self.displayMol.add_representation('ball+stick', selection='@1',
                                                       color='blue', opacity=0)  # Representation 2 = co_insertable
                    self.displayMol.add_representation('ball+stick', selection='@1',
                                                       color='green', opacity=0)  # Representation 3 = protrusion

                    self.displayMol.layout.width = '600px'
                    self.displayMol.layout.height = '600px'
                    # self.displayMol._set_selection('@', repr_index=1)
                    self.tab1.children = (list(self.tab1.children)[0], self.displayMol)

                    # Create a copy of the pdb so request are faster
                    self.parent.pdb = self.parent.dataset.query('cathpdb == @self.pdbWidget.value')
                    self.parent.haspdb = True
                    # Update exposition and entropy values

                    rangeExposition = (self.parent.pdb.sasa_rel_dssp.min(), self.parent.pdb.sasa_rel_dssp.max())
                    self.expositionWidget.min = rangeExposition[0]
                    self.expositionWidget.max = rangeExposition[1]
                    self.expositionWidget.value = rangeExposition[1]

                    rangeEntropy = (self.parent.pdb.shannonH10.min(), self.parent.pdb.shannonH10.max())
                    if np.isnan(rangeEntropy[0]) and np.isnan(rangeEntropy[1]):
                        self.entropyWidget.disable = True
                    else:
                        self.entropyWidget.disable = False
                        self.entropyWidget.min = rangeEntropy[0]
                        self.entropyWidget.max = rangeEntropy[1]
                        self.entropyWidget.value = rangeEntropy[1]

                    # self.out.clear_output()

        def set_entropy_as_bfactor(self, pdbname):
            entropy = self.parent.dataset.query("cathpdb == @pdbname").shannonH10.to_numpy()
            entropy = np.nan_to_num(entropy, nan=-1)
            ptp = np.ptp(
                entropy[~np.isnan(entropy)])  # copy of entropy without nan. ptp is the range that takes values.
            entropy = ((entropy - np.nanmin(entropy)) / ptp)
            self.mdaMol.atoms.tempfactors = entropy

        def update_selection(self, change=None):
            with self.out:
                # exposition
                exposed = self.parent.pdb.query(
                    f"sasa_rel_dssp {self.expositionSignWidget.value} @self.expositionWidget.value").residue_number.unique()
                # exposition
                conserved = self.parent.pdb.query(
                    f"shannonH10 {self.entropySignWidget.value} @self.entropyWidget.value").residue_number.unique()

                # Add segment
                IncludeResid_segments = []
                for segment in self.IncludeSegments:
                    low = segment.children[0].value
                    up = segment.children[1].value
                    if low != 0 and up != 0:
                        IncludeResid_segments.extend(list(range(low,
                                                                up)))
                # Exclude segment
                ExcludeResid_segments = []
                for segment in self.ExcludeSegments:
                    low = segment.children[0].value
                    up = segment.children[1].value
                    if low != 0 and up != 0:
                        ExcludeResid_segments.extend(list(range(low,
                                                                up)))
                # take common amino acids between 3 lists
                s1 = set(exposed)
                s2 = set(conserved)
                s3 = set(IncludeResid_segments)
                set1 = s1.intersection(s2)
                selected_resid = list(set1.intersection(s3))

                selection_string = ' or '.join(f"resid {x}" for x in selected_resid)
                if len(ExcludeResid_segments) > 0:
                    selection_string = '(' + selection_string + ')' + "and not (" + ' or '.join(
                        f"resid {x}" for x in ExcludeResid_segments) + ')'

                # Remove N,C,O
                # selection_string = "not type H and (" + selection_string + ")" #+" and not (name C or name O or name N)"
                ids = self.mdaMol.select_atoms(selection_string).ix
                self.displayMol._set_selection("not hydrogen and sidechainAttached and @" + ','.join(map(str, ids)),
                                               repr_index=1)
                self.displayMol.update_representation(component=0, repr_index=1, colorScheme='element')

        def on_createConvexhull_button_clicked(self, change):
            def normal(self, points, face):
                u = self.vect(points, face[0], face[1])
                v = self.vect(points, face[0], face[-1])
                return self.cross(u, v)

            with self.out:

                positions = self.mdaMol.select_atoms("name CB").positions
                hull = ConvexHull(positions)

                mesh = []
                for triangles in hull.simplices[hull.good]:
                    triangles.sort()
                    vertices = []
                    for vertex in triangles:
                        vertices.extend(positions[vertex].tolist())
                    mesh.append(vertices)
                mesh = np.asarray(mesh).flatten().tolist()
                color = [[0, 0, 1]] * len(mesh)  # RGB, for now let's fix it to blue
                color = np.asarray(color).flatten().tolist()
                shape = self.displayMol.shape
                # shape.add_mesh(mesh,color)

                shape.add_mesh(
                    mesh,
                    color,
                )

                self.displayMol.update_representation(component=1, repr_index=0, opacity=0.5, side="double", )

        def on_coinsertable_ticked(self, change):
            with self.out:
                if self.parent.haspdb:
                    if self.showCoinsertableWidget.value == True:
                        selected_resid = self.parent.pdb.query("is_co_insertable == 1").residue_number.unique()
                        selection_string = ' or '.join(f"resid {x}" for x in selected_resid)
                        ids = self.mdaMol.select_atoms(selection_string).ids
                        self.displayMol._set_selection(".CB and @" + ','.join(map(str, ids)), repr_index=2)
                        self.displayMol.update_representation(component=0, repr_index=2, opacity=1, radiusScale=5)
                    else:
                        self.displayMol._set_selection('@2', repr_index=2)
                        self.displayMol.update_representation(component=0, repr_index=2, opacity=0)

        def on_protrusion_ticked(self, change):
            with self.out:
                if self.parent.haspdb:
                    if self.showProtrusionWidget.value == True:
                        selected_resid = self.parent.pdb.query("is_hydrophobic_protrusion == 1").residue_number.unique()
                        selection_string = ' or '.join(f"resid {x}" for x in selected_resid)
                        ids = self.mdaMol.select_atoms(selection_string).ids
                        self.displayMol._set_selection(".CB and @" + ','.join(map(str, ids)), repr_index=3)
                        self.displayMol.update_representation(component=0, repr_index=3, opacity=1, radiusScale=5)

                    else:
                        self.displayMol._set_selection('@3', repr_index=3)
                        self.displayMol.update_representation(component=0, repr_index=3, opacity=0)

        def on_residueAnalysis_button_clicked(self, change):
            with self.out:
                self.parent.plot_residue_composition(self.parent.ibs, self.parent.nonibs)

        def on_run_button_clicked(self, change):
            with self.out:
                self.out.clear_output()

                domain = self.domainsWidget.value

                df = self.parent.dataset.query("domain == @domain")
                pdb = self.pdbWidget.value

                IncludeResidueSegments = []
                for segment in self.IncludeSegments:
                    start = segment.children[0].value
                    end = segment.children[1].value
                    if start != 0 and end != 0:
                        IncludeResidueSegments.append([start, end])

                # Exclude segment
                ExcludeResidueSegments = []
                for segment in self.ExcludeSegments:
                    start = segment.children[0].value
                    end = segment.children[1].value
                    if start != 0 and end != 0:
                        ExcludeResidueSegments.append([start, end])

                extendSS = self.extendSSWidget.value
                alignment = self.alignmentWidget.value
                onlyCoilSS = self.onlyCoilSSWidget.value
                cathCluster = self.cathClusterWidget.value
                uniref = self.unirefWidget.value
                excludeStrand = self.excludeStrand.value


                self.parent.tag_ibs(dataset=df,
                                    domain=domain,
                                    pdbreference=pdb,
                                    includeResidueRange=IncludeResidueSegments,
                                    excludeResidueRange=ExcludeResidueSegments,
                                    extendSS=extendSS,
                                    withAlignment=alignment,
                                    onlyC=onlyCoilSS,
                                    cathCluster=cathCluster,
                                    Uniref=uniref,
                                    addSequence=True,  # Always true for now,
                                    extendAlign=True,  # Always true for now,
                                    excludeStrand=excludeStrand,
                                    )

                # # Put dataset in the parent class
                # self.parent.ibs = tempout[0]
                # self.parent.nonibs = tempout[1]

                print("> IBS/NON-IBS datasets generated, you can run analysis now")
                self.residueAnalysisButton.disabled = False
            # with output:
            #    print(domain,pdb,residueRange,extendSS, alignment, onlyCoilSS, cathCluster, uniref)