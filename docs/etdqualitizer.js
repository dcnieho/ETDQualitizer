console.log("loading pyodide...");
async function main() {

  pyodide = await loadPyodide({
    indexURL: "https://cdn.jsdelivr.net/pyodide/v0.27.3/full/",
  });
  await pyodide.loadPackage("pandas");
  console.log("Pyodide/python loaded");
  document.getElementById('pyodideinfo').innerHTML = "";
  document.getElementById('pyodideinfo').style.display = 'none';

  await pyodide.loadPackage("micropip");
  const micropip = pyodide.pyimport("micropip");

  try {
    await micropip.install('ETDQualitizer');
    console.log("ETDQualitizer loaded");
    pyodide.runPython(`
        from importlib.metadata import version
        etdqver = version('ETDQualitizer')
    `);
    etdqver = pyodide.globals.get('etdqver');
    document.title = document.title + ' v' + etdqver;
    document.getElementById('etdqTitle').innerHTML = document.getElementById('etdqTitle').innerHTML + ' v' + etdqver;
    document.getElementById('etdqVersion').innerHTML = 'Version loaded: ' + etdqver;
    document.getElementById('bestand').disabled = false;
    document.getElementById('lblbestand').style.backgroundColor ='cadetblue';
    document.getElementById('lblbestand').style.color ='aliceblue';
    document.getElementById('lijst').innerHTML = document.getElementById('lijst').innerHTML + ' please select one or more files from your PyschoPy procedure, enter your<br/>setup geometry below and click Do it.';
    document.getElementById('lblbestand').style.backgroundColor ='cadetblue';
  } catch (error) {
    document.getElementById('pyodideinfo').innerHTML = document.getElementById('pyodideinfo').innerHTML +  " and ETDQualitizer failed to install";
    console.error(error);
  }

}

main();
