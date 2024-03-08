function getSystemInfo() {

  let half1 = api.getElements([...Array(elements.size/2).keys()]);
  let half2 = api.getElements(Array.from({length: elements.size - half1.length}, (_, i) => i + half1.length));
  colorElements(new THREE.Color(0x0000ff), half1)

  // label first half of elements
  colorElements(new THREE.Color(0x8effff), half2)
  let half1_s = new Set(half1)
  let half2_s = new Set(half2)

  // identify pairing information.  I'm using oxView here rather than output bonds because its more permissive of not-quite-bonded nucleotides and so is less likley to miss nucleotides that are randomly frayed in the particular simulation step I used for identifying these areas.
  findBasepairs(3)

  // identify the pairs of our second set (remembering that oxDNA is 3' -> 5')
  let pairs1 = []
  half2.forEach(e => { if (e.pair) { pairs1.push(e.pair) } })
  pairs1_s = new Set(pairs1)

  //set intersection function from MDN docs
  function intersection(setA, setB) { 
      let _intersection = new Set() 
      for (let elem of setB) { 
          if (setA.has(elem)) { 
              _intersection.add(elem) 
          } 
      } 
      return _intersection 
  }

  // identify points where half1 binds with itself
  let badNucs = intersection(pairs1_s, half2_s);
  badNucs = Array.from(badNucs);

  colorElements(new THREE.Color(0xff0000), badNucs);

  return [half2, badNucs];
}

///////////////////////////////////////////////////////////////////////////////////////////
// Run getSystemInfo()
///////////////////////////////////////////////////////////////////////////////////////////

half2_badNucs = getSystemInfo()
half2 = half2_badNucs[0]
badNucs = half2_badNucs[1]

///////////////////////////////////////////////////////////////////////////////////////////
// At this point we now modify the 3` and 5` ends of the template by hand in oxView
///////////////////////////////////////////////////////////////////////////////////////////

// Data to modify template by hand
unstruct_3 = 'GTGGTGGGGCCACTTGTCGAGGAGCGGGAACGAGTGGTA' // first select 3` end and extend
utr_5 = 'ccaccgcccaagagagactcagacacccctggtcttcttatgatcaaataa' // then select 3` end and extend the 5` 
utr_3 = 'tgagctggagcctcggtggcctagcttcttgccccttgggcctccccccagcccctcctccccttcctgcacccgtacccccgtggtctttgaataaagtctgagtgggcggc'
poly_a = 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'

/////////////////////////////////////////////////////////////////////////
// done with hand modification
/////////////////////////////////////////////////////////////////////////

function editSequence(half2, badNucs) {

  // The coding sequence of our gene in 3` to 5` that will be structred
  let seq = 'CCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAATGATAATAGACCGGT'

  //we previously selected them 3' -> 5', so we need the reverse.
  half2 = half2.reverse()

  //shoulda just set this the first time
  let badNucs_s = new Set(badNucs)

  //run through the list and highlight all the nucleotides we've changed.
  half2.forEach((e, i) => {
    if (i < seq.length) {
      e.setType(seq[i])
      if (!badNucs_s.has(e) && e.pair) {
        e.pair.setType(e.getComplementaryType())
      }
      api.selectElements([e], true)
    }
  })

  //force the colors to update because I can't remember the instance call for colors.
  api.toggleBaseColors()
  
  return null;
}

///////////////////////
// run this and thecn change the color of the selected bases in oxView

editSequence(half2, badNucs)

//////////////////////////////////////////////////////////////////////////////
// At this point we select the red nucleotides by color in oxView and run the following code
//////////////////////////////////////////////////////////////////////////////
selectedBases.forEach((e) => {
    e.setType(e.pair.getComplementaryType())
  
})

api.toggleBaseColors()