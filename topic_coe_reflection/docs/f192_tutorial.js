/* 纯本地阅读辅助：目录定位、搜索、字号、复制及图片放大。 */
(() => {
  'use strict';
  document.body.classList.add('js-ready');
  const chapters = [...document.querySelectorAll('main > header, main > section')];
  const nav = document.querySelector('.book-nav');
  const navLinks = [...nav.querySelectorAll('a[href^="#"]')];
  const location = document.querySelector('.current-location');
  const toc = document.querySelector('.chapter-toc');
  const progress = document.querySelector('.progress-track span');
  const menu = document.querySelector('.menu-button');
  menu.addEventListener('click', () => {
    const open = document.body.classList.toggle('nav-open');
    menu.setAttribute('aria-expanded', String(open));
  });
  navLinks.forEach(a => a.addEventListener('click', () => {
    document.body.classList.remove('nav-open');
    menu.setAttribute('aria-expanded', 'false');
  }));
  let bodySize = 17;
  const smaller = document.querySelector('[data-size="smaller"]');
  const larger = document.querySelector('[data-size="larger"]');
  function setSize(delta) {
    bodySize = Math.max(15, Math.min(21, bodySize + delta));
    document.documentElement.style.setProperty('--body-size', `${bodySize}px`);
    smaller.disabled = bodySize === 15;
    larger.disabled = bodySize === 21;
    schedule();
  }
  smaller.addEventListener('click', () => setSize(-1));
  larger.addEventListener('click', () => setSize(1));
  const search = document.querySelector('#chapter-search');
  const results = document.querySelector('.search-results');
  const index = chapters.map(c => ({id:c.id, title:c.querySelector('h1,h2').textContent.replace(/\s+/g,' ').trim(), text:c.textContent.toLowerCase()}));
  search.addEventListener('input', () => {
    results.replaceChildren();
    const query = search.value.trim().toLowerCase();
    if (!query) return;
    const hits = index.filter(c => c.text.includes(query));
    const label = document.createElement('div');
    label.textContent = hits.length ? `找到 ${hits.length} 个章节` : '没有匹配章节，试试“曝光”或“reshape”';
    results.append(label);
    hits.forEach(hit => {
      const a = document.createElement('a'); a.href = `#${hit.id}`; a.textContent = hit.title;
      a.addEventListener('click', () => {document.body.classList.remove('nav-open');menu.setAttribute('aria-expanded','false');});
      results.append(a);
    });
  });
  // 生成行号时保留原始字符串；复制时不用视觉行号或HTML。
  document.querySelectorAll('pre > code').forEach((code, index) => {
    const source = code.textContent;
    const toolbar = document.createElement('div'); toolbar.className = 'code-toolbar';
    const label = document.createElement('span');label.textContent = `MATLAB · 代码 ${String(index + 1).padStart(2,'0')}`;
    const button = document.createElement('button');button.type='button';button.textContent='复制代码';
    button.setAttribute('aria-label',`复制第 ${index+1} 个 MATLAB 代码块`);
    toolbar.append(label,button);code.parentElement.before(toolbar);
    code.replaceChildren();code.classList.add('numbered-code');
    source.split('\n').forEach((line,i,all) => {
      const span = document.createElement('span');span.className='code-line';span.dataset.line=String(i+1);
      if (line.trimStart().startsWith('%')) span.classList.add('code-comment');
      span.textContent=line;code.append(span);
      if(i<all.length-1)code.append(document.createTextNode('\n'));
    });
    button.addEventListener('click', async () => {
      try {
        await navigator.clipboard.writeText(source);button.textContent='已复制';
      } catch (_) {
        const selection=window.getSelection(),range=document.createRange();
        range.selectNodeContents(code);selection.removeAllRanges();selection.addRange(range);
        button.textContent='已选中，请按 Ctrl+C';
      }
      setTimeout(()=>{button.textContent='复制代码';},2500);
    });
  });
  document.querySelectorAll('figure img').forEach(img => {
    const actions=document.createElement('div');actions.className='figure-actions';
    const button=document.createElement('button');button.type='button';button.textContent='查看大图';
    actions.append(button);img.after(actions);
    button.addEventListener('click',()=>{
      const dialog=document.createElement('dialog');dialog.className='image-dialog';
      const top=document.createElement('div'),title=document.createElement('span'),close=document.createElement('button');
      title.textContent=img.alt;close.type='button';close.textContent='关闭';
      const large=document.createElement('img');large.src=img.currentSrc||img.src;large.alt=img.alt;
      top.append(title,close);dialog.append(top,large);document.body.append(dialog);
      close.addEventListener('click',()=>dialog.close());
      dialog.addEventListener('click',e=>{if(e.target===dialog)dialog.close();});
      dialog.addEventListener('close',()=>{dialog.remove();button.focus();});
      dialog.showModal();
    });
  });
  let current = null, scheduled = false;
  function update() {
    scheduled=false;
    let active=chapters[0];
    for(const c of chapters){if(c.getBoundingClientRect().top<=140)active=c;else break;}
    const isEnd=window.scrollY+window.innerHeight>=document.documentElement.scrollHeight-8;
    if(isEnd)active=chapters[chapters.length-1];
    if(current!==active){
      current=active;
      navLinks.forEach(a=>{if(a.hash===`#${active.id}`)a.setAttribute('aria-current','location');else a.removeAttribute('aria-current');});
      const activeLink=navLinks.find(a=>a.hash===`#${active.id}`);
      if(activeLink&&getComputedStyle(nav).position==='fixed'){
        const rect=activeLink.getBoundingClientRect(),bounds=nav.getBoundingClientRect();
        if(rect.top<bounds.top+30||rect.bottom>bounds.bottom-30)nav.scrollTop+=rect.top-bounds.top-120;
      }
      const heading=active.querySelector('h1,h2').textContent.replace(/\s+/g,' ').trim();
      location.textContent=heading;
      toc.replaceChildren();
      const title=document.createElement('h2');title.textContent='本章导航';
      const name=document.createElement('p');name.textContent=heading;toc.append(title,name);
      active.querySelectorAll('h3').forEach(h=>{const a=document.createElement('a');a.href=`#${h.id}`;a.textContent=h.textContent;toc.append(a);});
    }
    const heads=[...active.querySelectorAll('h3')];
    let selected=heads[0];heads.forEach(h=>{if(h.getBoundingClientRect().top<=155)selected=h;});
    toc.querySelectorAll('a').forEach(a=>{if(selected&&a.hash===`#${selected.id}`)a.setAttribute('aria-current','location');else a.removeAttribute('aria-current');});
    const height=document.documentElement.scrollHeight-window.innerHeight;
    progress.style.width=`${height>0?Math.max(0,Math.min(100,window.scrollY/height*100)):0}%`;
  }
  function schedule(){if(!scheduled){scheduled=true;window.requestAnimationFrame(update);}}
  window.addEventListener('scroll',schedule,{passive:true});
  window.addEventListener('resize',schedule);
  document.querySelectorAll('img').forEach(img=>img.addEventListener('load',schedule));
  // 打印时展开答案，打印结束恢复用户此前的折叠状态。
  let printState=[];
  window.addEventListener('beforeprint',()=>{printState=[...document.querySelectorAll('details')].map(d=>[d,d.open]);printState.forEach(([d])=>{d.open=true;});});
  window.addEventListener('afterprint',()=>{printState.forEach(([d,open])=>{d.open=open;});});
  update();
})();
