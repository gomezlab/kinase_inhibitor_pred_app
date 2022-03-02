function button_reaction() {
	in_process.style.display = 'block'; 
	instructions.style.display = 'hide';
	document.getElementById("submit_random_geo").disabled = true;
	document.getElementById("submit_geo").disabled = true;
	document.getElementById("submit_upload_seq").disabled = true;
}