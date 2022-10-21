
Network construction methods
============================

.. raw:: html

	<table border="1" class="dataframe">
	  <thead>
	    <tr>
	      <th>Residue information extracted from trajectories</th>
	      <th>Package</th>
	      <th>Correlation measurement</th>
	      <th>Atom/angle tracked</th>
	      <th></th>
	    </tr>
	  </thead>
	  <tbody>
	    <tr>
	      <th rowspan="8">Atoms' movement correlation</th>
	      <th>MD-TASK</th>
	      <th>Pearson's</th>
	      <th>Carbon α</th>
	      <td>MDTASK</td>
	    </tr>
	    <tr>
	      <th rowspan="2">pytraj</th>
	      <th rowspan="2">Pearson's</th>
	      <th>Carbon α</th>
	      <td>pytraj_CA</td>
	    </tr>
	    <tr>
	      <th>Carbon β</th>
	      <td>pytraj_CB</td>
	    </tr>
	    <tr>
	      <th>dynetan</th>
	      <th>MI</th>
	      <th>Whole residue</th>
	      <td>dynetan</td>
	    </tr>
	    <tr>
	      <th rowspan="4">correlationplus</th>
	      <th rowspan="2">Pearson's</th>
	      <th>Carbon α</th>
	      <td>correlationplus_CA_Pear</td>
	    </tr>
	    <tr>
	      <th>Residue COM</th>
	      <td>correlationplus_COM_Pear</td>
	    </tr>
	    <tr>
	      <th rowspan="2">LMI</th>
	      <th>Carbon α</th>
	      <td>correlationplus_CA_LMI</td>
	    </tr>
	    <tr>
	      <th>Residue COM</th>
	      <td>correlationplus_COM_LMI</td>
	    </tr>
	    <tr>
	      <th rowspan="68">Dihedrals' movement correlation</th>
	      <th rowspan="4">correlationplus</th>
	      <th rowspan="4">Pearson's</th>
	      <th>Phi</th>
	      <td>correlationplus_Phi</td>
	    </tr>
	    <tr>
	      <th>Psi</th>
	      <td>correlationplus_Psi</td>
	    </tr>
	    <tr>
	      <th>All backbone dihedrals (Phi and psi) (average)</th>
	      <td>correlationplus_Backbone_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>All backbone dihedrals (Phi and psi) (max. value)</th>
	      <td>correlationplus_Backbone_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th rowspan="12">AlloViz</th>
	      <th rowspan="12">MI</th>
	      <th>Phi</th>
	      <td>AlloViz_Phi</td>
	    </tr>
	    <tr>
	      <th>Psi</th>
	      <td>AlloViz_Psi</td>
	    </tr>
	    <tr>
	      <th>All backbone dihedrals (Phi and psi) (average)</th>
	      <td>AlloViz_Backbone_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>All backbone dihedrals (Phi and psi) (max. value)</th>
	      <td>AlloViz_Backbone_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>Chi1</th>
	      <td>AlloViz_Chi1</td>
	    </tr>
	    <tr>
	      <th>Chi2</th>
	      <td>AlloViz_Chi2</td>
	    </tr>
	    <tr>
	      <th>Chi3</th>
	      <td>AlloViz_Chi3</td>
	    </tr>
	    <tr>
	      <th>Chi4</th>
	      <td>AlloViz_Chi4</td>
	    </tr>
	    <tr>
	      <th>All side-chain dihedrals (average)</th>
	      <td>AlloViz_Sidechain_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>All side-chain dihedrals (max. value)</th>
	      <td>AlloViz_Sidechain_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>All dihedrals (average)</th>
	      <td>AlloViz_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>All dihedrals (max. value)</th>
	      <td>AlloViz_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th rowspan="48">CARDS</th>
	      <th>MI</th>
	      <th>Phi</th>
	      <td>CARDS_MI_Phi</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>Phi</th>
	      <td>CARDS_Disorder_Phi</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>Phi</th>
	      <td>CARDS_Disorder_mediated_Phi</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>Phi</th>
	      <td>CARDS_Holistic_Phi</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>Psi</th>
	      <td>CARDS_MI_Psi</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>Psi</th>
	      <td>CARDS_Disorder_Psi</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>Psi</th>
	      <td>CARDS_Disorder_mediated_Psi</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>Psi</th>
	      <td>CARDS_Holistic_Psi</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>All backbone dihedrals (Phi and psi) (average)</th>
	      <td>CARDS_MI_Backbone_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>All backbone dihedrals (Phi and psi) (average)</th>
	      <td>CARDS_Disorder_Backbone_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>All backbone dihedrals (Phi and psi) (average)</th>
	      <td>CARDS_Disorder_mediated_Backbone_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>All backbone dihedrals (Phi and psi) (average)</th>
	      <td>CARDS_Holistic_Backbone_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>All backbone dihedrals (Phi and psi) (max. value)</th>
	      <td>CARDS_MI_Backbone_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>All backbone dihedrals (Phi and psi) (max. value)</th>
	      <td>CARDS_Disorder_Backbone_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>All backbone dihedrals (Phi and psi) (max. value)</th>
	      <td>CARDS_Disorder_mediated_Backbone_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>All backbone dihedrals (Phi and psi) (max. value)</th>
	      <td>CARDS_Holistic_Backbone_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>Chi1</th>
	      <td>CARDS_MI_Chi1</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>Chi1</th>
	      <td>CARDS_Disorder_Chi1</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>Chi1</th>
	      <td>CARDS_Disorder_mediated_Chi1</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>Chi1</th>
	      <td>CARDS_Holistic_Chi1</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>Chi2</th>
	      <td>CARDS_MI_Chi2</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>Chi2</th>
	      <td>CARDS_Disorder_Chi2</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>Chi2</th>
	      <td>CARDS_Disorder_mediated_Chi2</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>Chi2</th>
	      <td>CARDS_Holistic_Chi2</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>Chi3</th>
	      <td>CARDS_MI_Chi3</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>Chi3</th>
	      <td>CARDS_Disorder_Chi3</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>Chi3</th>
	      <td>CARDS_Disorder_mediated_Chi3</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>Chi3</th>
	      <td>CARDS_Holistic_Chi3</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>Chi4</th>
	      <td>CARDS_MI_Chi4</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>Chi4</th>
	      <td>CARDS_Disorder_Chi4</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>Chi4</th>
	      <td>CARDS_Disorder_mediated_Chi4</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>Chi4</th>
	      <td>CARDS_Holistic_Chi4</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>All side-chain dihedrals (average)</th>
	      <td>CARDS_MI_Sidechain_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>All side-chain dihedrals (average)</th>
	      <td>CARDS_Disorder_Sidechain_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>All side-chain dihedrals (average)</th>
	      <td>CARDS_Disorder_mediated_Sidechain_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>All side-chain dihedrals (average)</th>
	      <td>CARDS_Holistic_Sidechain_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>All side-chain dihedrals (max. value)</th>
	      <td>CARDS_MI_Sidechain_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>All side-chain dihedrals (max. value)</th>
	      <td>CARDS_Disorder_Sidechain_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>All side-chain dihedrals (max. value)</th>
	      <td>CARDS_Disorder_mediated_Sidechain_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>All side-chain dihedrals (max. value)</th>
	      <td>CARDS_Holistic_Sidechain_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>All dihedrals (average)</th>
	      <td>CARDS_MI_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>All dihedrals (average)</th>
	      <td>CARDS_Disorder_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>All dihedrals (average)</th>
	      <td>CARDS_Disorder_mediated_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>All dihedrals (average)</th>
	      <td>CARDS_Holistic_Dihs_Avg</td>
	    </tr>
	    <tr>
	      <th>MI</th>
	      <th>All dihedrals (max. value)</th>
	      <td>CARDS_MI_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>Pure-disorder MI</th>
	      <th>All dihedrals (max. value)</th>
	      <td>CARDS_Disorder_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>Disorder-mediated MI</th>
	      <th>All dihedrals (max. value)</th>
	      <td>CARDS_Disorder_mediated_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th>Holistic MI</th>
	      <th>All dihedrals (max. value)</th>
	      <td>CARDS_Holistic_Dihs_Max</td>
	    </tr>
	    <tr>
	      <th rowspan="4">MDEntropy</th>
	      <th rowspan="4">MI</th>
	      <th>Phi</th>
	      <td>MDEntropy_Phi</td>
	    </tr>
	    <tr>
	      <th>Psi</th>
	      <td>MDEntropy_Psi</td>
	    </tr>
	    <tr>
	      <th>All backbone dihedrals (Phi and psi) (average)</th>
	      <td>MDEntropy_Dihs</td>
	    </tr>
	    <tr>
	      <th>Alpha angle</th>
	      <td>MDEntropy_AlphaAngle</td>
	    </tr>
	    <tr>
	      <th rowspan="6">Contacts</th>
	      <th>MDEntropy</th>
	      <th>MI</th>
	      <th>Contact frequency</th>
	      <td>MDEntropy_Contacts</td>
	    </tr>
	    <tr>
	      <th>GetContacts</th>
	      <th>None</th>
	      <th>Contact frequency</th>
	      <td>GetContacts</td>
	    </tr>
	    <tr>
	      <th rowspan="3">PyInteraph2</th>
	      <th rowspan="3">None</th>
	      <th>Contact frequency</th>
	      <td>PyInteraph2_Atomic_Contacts_Occurrence</td>
	    </tr>
	    <tr>
	      <th>Contact strength</th>
	      <td>PyInteraph2_Atomic_Contacts_Strength</td>
	    </tr>
	    <tr>
	      <th>Residue COM contacts</th>
	      <td>PyInteraph2_COM_Contacts</td>
	    </tr>
	    <tr>
	      <th>PyInteraph2 (with Rg correction)</th>
	      <th>None</th>
	      <th>Residue COM contacts</th>
	      <td>PyInteraph2_COM_Contacts_Corrected</td>
	    </tr>
	    <tr>
	      <th rowspan="3">Interaction energy</th>
	      <th>PyInteraph2</th>
	      <th>None</th>
	      <th>Whole residue</th>
	      <td>PyInteraph2_Energy</td>
	    </tr>
	    <tr>
	      <th rowspan="2">gRINN</th>
	      <th>None</th>
	      <th>Whole residue</th>
	      <td>gRINN</td>
	    </tr>
	    <tr>
	      <th>Pearson's</th>
	      <th>Whole residue</th>
	      <td>gRINN_corr</td>
	    </tr>
	  </tbody>
	</table>
